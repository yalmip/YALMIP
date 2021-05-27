function output = bmibnb(p)
%BMIBNB Global solver based on branch-and-bound
%
% BMIBNB is never called by the user directly, but is called by YALMIP from
% OPTIMIZE, by choosing the solver tag 'bmibnb' in sdpsettings
%
% The behaviour of BMIBNB can be altered using the fields in the field
% 'bmibnb' in SDPSETTINGS
%
% bmibnb.lowersolver    - Solver for lower bound [solver tag ('')]
% bmibnb.uppersolver    - Solver for upper bound [solver tag ('')]
% bmibnb.lpsolver       - Solver for LP bound tightening [solver tag ('')]
%
% bmibnb.maxiter        - Maximum number nodes [int (100)]
% bmibnb.maxtime        - Maximum CPU time (sec.) [int (3600)]
%
% bmibnb.relgaptol      - Tolerance on relative objective (in %) error (UPPER-LOWER)/(1+|UPPER|) [real (0.01)]
% bmibnb.absgaptol      - Tolerance on objective error (UPPER-LOWER) [real (0.01)]
%
% bmibnb.target         - Exit if upper bound<=target [double (-inf)]
% bmibnb.lowertarget    - Exit if lower bound>=target [double ( inf)]
%
% bmibnb.branchmethod   - Branch strategy ['maxvol' | 'best' ('best')]
% bmibnb.branchrule     - Branch position ['omega' | 'bisect' ('omega')]
%
% bmibnb.lpreduce       - Improve variable bounds using LP (-1 BMIBNB decides, 0 means no reduction, 1 means all variables)
% bmibnb.lowrank        - Partition variables into two disjoint sets and branch on smallest [ 0|1 (0)]
% bmibnb.roottight      - Improve variable bounds in root using full problem (-1 BMIBNB decides, 0 off, 1 on)
% bmibnb.vartol         - Cut tree when x_U-x_L < vartol on all branching variables
% bmibnb.pdtol          - A number is declared non-negative if larger than...[ double (1e-6)]
% bmibnb.eqtol          - A number is declared zero if abs(x) smaller than...[ double (1e-6)]

showprogress('Branch and bound started',p.options.showprogress);

% *************************************************************************
% INITIALIZE DIAGNOSTICS AND DISPLAY LOGICS (LOWER VERBOSE IN SOLVERS)
% *************************************************************************
switch max(min(p.options.verbose,3),0)
case 0
    p.options.bmibnb.verbose = 0;
case 1
    p.options.bmibnb.verbose = 1;
    p.options.verbose = 0;
case 2
    p.options.bmibnb.verbose = 2;
    p.options.verbose = 0;
case 3
    p.options.bmibnb.verbose = 2;
    p.options.verbose = 1;
otherwise
    p.options.bmibnb.verbose = 0;
    p.options.verbose = 0;
end

if p.options.bmibnb.verbose>0
    disp('* Starting YALMIP global branch & bound.');
    disp(['* Upper solver     : ' p.solver.uppersolver.tag]);
    disp(['* Lower solver     : ' p.solver.lowersolver.tag]);    
    disp(['* LP solver        : ' p.solver.lpsolver.tag]); 
end

% CPLEX handles options insanely slow, so we remove all default values in
% solver options
p.options = pruneOptions(p.options);

% We will fake QP-support in Mosek SDP
if strcmpi(p.solver.lowersolver.tag,'mosek')
    if any(p.K.s)
       p.solver.lowersolver.objective.quadratic.convex = 1;
    end
end

timing.total = tic;
timing.uppersolve = 0;
timing.lowersolve = 0;
timing.lpsolve = 0;
timing.heuristics = 0;
if ~isempty(p.F_struc)
    if any(isnan(p.F_struc) | isinf(p.F_struc))
        output = createOutputStructure(1);        
        output.solved_nodes = 0;
        output.Primal       = zeros(length(p.c),1);
        output.infostr      = yalmiperror(output.problem,'BNB');
        output.solvertime   = 0;
        output.lower = inf;
        return
    end
end

% Make the evaluation structure a bit more unified
% by embedding all trailing arguments and creating
% short-hand propert to evaluate function value
for i = 1:length(p.evalMap)
    p.evalMap{i}.properties.function = @(x)(feval(p.evalMap{i}.fcn,x,p.evalMap{i}.arg{2:end-1}));
end

p.nonshiftedQP.Q =[];
p.nonshiftedQP.c =[];
p.nonshiftedQP.f =[];

% *************************************************
% Decide on strategy for nonlinear SDP cone
% *************************************************
if p.options.bmibnb.uppersdprelax < 0
    if p.solver.uppersolver.constraint.inequalities.semidefinite.polynomial
        % No reason to add nonlinear cuts as we have a nonlinear SDP solver
        p.options.bmibnb.uppersdprelax = 0;
    else
        %Solver does not support SDP cone, so accept iterations
        p.options.bmibnb.uppersdprelax = -p.options.bmibnb.uppersdprelax;
    end
elseif p.options.bmibnb.uppersdprelax > 0
    % Remove knowledge of SP cone from upper solver
    p.solver.uppersolver.constraint.inequalities.semidefinite.polynomial = 0;
    p.solver.uppersolver.constraint.inequalities.semidefinite.quadratic = 0;
    p.solver.uppersolver.constraint.inequalities.semidefinite.linear = 0;    
end

% *************************************************************************
% Assume feasible (this property can be changed in the presolve codes
% *************************************************************************
p.feasible = 1;

% *************************************************************************
% Ensure that bound propagation is performed
% *************************************************************************
p.changedbounds = 1;

% *************************************************************************
% Initialize some diagnostic counters
% *************************************************************************
p.counter.lpsolved = 0;
p.counter.lowersolved = 0;
p.counter.uppersolved = 0;
p.counter.heuristics = 0;

% *************************************************************************
% Extend partially specified initials
% *************************************************************************
p = completeInitial(p);

% *************************************************************************
% Some preprocessing to put the model in a better format
% *************************************************************************
p = convert_perspective_log(p);
n_in = length(p.c);

% *************************************************************************
% Copy positivity terms in cones to LPs to use these in bound propagation
% *************************************************************************
if nnz(p.K.q)>0
    pos = p.K.f + p.K.l + cumsum([1 p.K.q(1:end-1)]);
    p.F_struc = [p.F_struc(1:p.K.f+p.K.l,:);p.F_struc(pos,:);p.F_struc(1+p.K.f+p.K.l:end,:)];
    p.K.l = p.K.l + length(pos);
end

% *************************************************************************
% Add diagonal cuts from if we are going to use a cutting strategi for
% semidefinite constraints. No reason to spend nodes on finding these and
% we do it early to enable detection of variable bounds
% *************************************************************************
p = bounds_from_cones_to_lp(p);

% *************************************************************************
% Save information about the applicability of some bound propagation
% *************************************************************************
p.boundpropagation.sepquad = 1;

% *************************************************************************
% Extract bounds from model using direct information available
% *************************************************************************
bounds = yalmip('getbounds',p.used_variables);
if isequal(size(bounds),size(p.lb))
    p.lb = max(p.lb,bounds(:,1));
    p.ub = min(p.ub,bounds(:,2));
end
p = compile_nonlinear_table(p);

if p.options.bmibnb.verbose>0
	disp('* -Extracting bounds from model');   
end
p = detect_quadratic_disjoints(p);
p = detect_quadratic_disjoints_generalized(p);
p = detect_hiddendelayedconvex_sdp(p);
p = presolve_bounds_from_domains(p);
p = presolve_bounds_from_modelbounds(p,1);

% *************************************************************************
% Start some bound propagation
% *************************************************************************
if p.options.bmibnb.verbose>0
	disp('* -Perfoming root-node bound propagation');   
end
p = propagate_bounds_from_convex_quadratic_ball(p);
p = propagate_bounds_from_evaluations(p);
p = update_monomial_bounds(p);
p = presolve_sortrows(p);
p = propagate_bounds_from_equalities(p); 
p = propagate_bounds_from_evaluations(p);
p = propagate_bounds_from_separable_quadratic_equality(p);
p = update_monomial_bounds(p);
p = presolve_bounds_from_inequalities(p);
p = propagate_bounds_from_evaluations(p);
p = update_monomial_bounds(p);
p = propagate_bounds_from_convex_quadratic_ball(p);

% *************************************************************************
% For quadratic nonconvex programming over linear constraints, we
% diagonalize the problem to obtain less number of bilinear terms. Not
% theoretically any reason to do so, but practical performance is better
% If we are using diagonalize=-1, a choice will be made by comparing
% two lowersolves, and those also generate feasible solutions
% we we can use for initial upper bound
% *************************************************************************
x_min = zeros(length(p.c),1);
upper = inf;
[p,timing,x_min,upper] = diagonalize_quadratic_program(p,timing,x_min,upper);
if p.diagonalized
    p = compile_nonlinear_table(p);
end

% *************************************************************************
% Try to generate a feasible solution, by using avialable x0 (if usex0=1),
% or by trying the zero solution, or my trying the (lb+ub)/2 solution. Note
% that we do this here before the problem possibly is bilinearized, thus
% avoiding to introduce possibly complicating bilinear constraints
% *************************************************************************
[p,x_min,upper] = initializesolution(p,x_min,upper);
solution_hist = [];
if ~isinf(upper)   
    solution_hist = [solution_hist x_min(p.linears)];
    if p.options.bmibnb.verbose 
        disp(['* -Trivial solution constructed (objective ' num2str(upper) ')']);
    end
end

% *************************************************************************
% If the upper bound solver is capable of solving the original problem,
% without bilinearizing it first, we might try to get a local solution,
% possibly based on our initial crude solution
% *************************************************************************
if solver_can_solve(p.solver.uppersolver,p)
    p.high_monom_model = [];
    p = build_recursive_scheme(p);
    p = compile_nonlinear_table(p);
    p = preprocess_bilinear_bounds(p);
    p = propagate_bounds_from_evaluations(p);
    p = propagate_bounds_from_upper(p,upper);
    p = propagate_bounds_from_complementary(p);
    p = propagate_bounds_from_monomials(p);
    p = propagate_bounds_from_equalities(p);
    p = propagate_bounds_from_monomials(p);   
    p = propagate_bounds_from_monomials(p);    
    p = propagate_bounds_from_evaluations(p);    
    if p.options.bmibnb.verbose>0
        fprintf('* -Calling upper solver ');   
    end
    % Note that upper solver can add cuts to model if it is an SDP
    [upper_,x_min_,info_text,numglobals,timing,p] = solve_upper_in_node(p,p,x_min,upper,x_min,p.solver.uppersolver.call,'',0,timing,p.options.bmibnb.uppersdprelax);
    if upper_ < upper
        solution_hist = [solution_hist x_min_(p.linears)];
        upper = upper_;
        x_min = x_min_;              
        p = propagate_bounds_from_upper(p,upper);
        p = update_monomial_bounds(p);
        p = propagate_bounds_from_evaluations(p);
        p = propagate_bounds_from_monomials(p);
        p = presolve_bounds_from_inequalities(p);
        p = propagate_bounds_from_equalities(p);
        p = update_monomial_bounds(p);
        if p.options.bmibnb.verbose>0
            disp(['(found a solution, objective ' num2str(upper) ')']);            
        end        
    else
        if ~isinf(upper) && p.options.bmibnb.verbose>0
            disp('(no better solution found)');
        elseif p.options.bmibnb.verbose>0
            disp('(no solution found)');
        end
    end
    p.counter.uppersolved = p.counter.uppersolved + 1;
end
if isempty(p.x0)
    p.x0 = zeros(length(p.c),1);
end

% *************************************************************************
% Adaptively turn off LP-based propagation
% *************************************************************************
if p.options.bmibnb.lpreduce == -1
    if size(p.F_struc,1)==0 && isempty(p.evalMap) && all(p.variabletype <= 2) && isempty(p.binary_variables) && isempty(p.integer_variables) 
        % No constraints and no operators which might introduce any
        % interesting cuts, so LP-based propagation will only be driven by
        % bounds and quadratic objective, which never can improve anything
        % over the simply bound propagators
        p.options.bmibnb.lpreduce = 0;
    else
        p.options.bmibnb.lpreduce = 1;
    end
end
if p.options.bmibnb.roottight == -1
    if size(p.F_struc,1)==0 && isempty(p.evalMap) && all(p.variabletype <= 2) &&  isempty(p.binary_variables) && isempty(p.integer_variables) 
        % No constraints and no operators which might introduce any
        % interesting cuts, so LP-based propagation will only be driven by
        % bounds and quadratic objective, which never can improve anything
        % over the simply bound propagators
        p.options.bmibnb.roottight = 0;
    else
        p.options.bmibnb.roottight = 1;
    end
end
% *************************************************************************
% When saving history of solutions, we only want original linear terms
% *************************************************************************
p.originallinears = p.linears;

% *************************************************************************
% Sigmonial terms are converted to evaluation based expressions.
% *************************************************************************
p = convert_sigmonial_to_sdpfun(p);

% *************************************************************************
% The bilinear solver does not support non-quadratic models
% Hence, we need to bilinearize the problem. However, the upper bound
% solver might support general problems, so we should keep a copy of the
% original problem also in that case
% *************************************************************************
original_variables = find(p.variabletype == 0);
[p_bilin,changed] = convert_polynomial_to_quadratic(p);
if (p.solver.uppersolver.constraint.equalities.polynomial &  p.solver.uppersolver.objective.polynomial)
    p_bilin.originalModel = p;
    p = p_bilin;
else
    p = p_bilin;
    p.originalModel = p_bilin;
end
if changed
    p = compile_nonlinear_table(p);
end

p.EqualityConstraintState = ones(p.K.f,1);
p.InequalityConstraintState = ones(p.K.l,1);
p = propagate_bounds_from_arbitrary_quadratics(p);
    

% *************************************************************************
% Build an evaluation tree for computing nonlinear terms when running the
% upper-bound solver, typically fmincon. Operators have to be applied in
% a particular order, to cope with terms such as sin(cos(x^2)^2)
p.originalModel = build_recursive_scheme(p.originalModel);

% *************************************************************************
% Fast derivbation of upper and lower bounds also requires an evaluation
% order
p = build_recursive_scheme(p);

% *************************************************************************
% Pre-calc lists of linear/bilinear/nonlinear variables (we have bilineared
% the model now, so the old precompiled table could be wrong
% *************************************************************************
p = compile_nonlinear_table(p);

% *************************************************************************
% Select branch variables. We should branch on all variables that are
% involved in bilinear terms.
% *************************************************************************
p.branch_variables = decide_branch_variables(p);
p.branch_variables = setdiff(p.branch_variables,p.evalVariables);
p.branch_variables = intersect(p.branch_variables,original_variables);

if p.diagonalized
    if p.solver.lowersolver.objective.quadratic.convex
        Q =  compileQuadratic(p.c,p,2);
        convex = find(diag(Q)>=0);
        p.branch_variables = setdiff(p.branch_variables,convex);
    end
end
        
% *************************************************************************
% Tighten bounds (might be useful after bilinearization?)
% *************************************************************************
p = presolve_implied_integer(p);
p = preprocess_bilinear_bounds(p);
p = propagate_bounds_from_evaluations(p);
p = propagate_bounds_from_equalities(p);

% *************************************************************************
% Now reduce the branch variables by removing bilinear terms that only have
% been introduced due to a convex quadratic objective
% *************************************************************************
p = reduce_bilinear_branching_variables(p);

if p.options.bmibnb.verbose>0
    disp(['* -Branch-variables : ' num2str(length(p.branch_variables))]);  
end
% *************************************************************************
% Simple pre-solve loop. The loop is needed for variables defined as w =
% x*y, x = t*u,y=..db
% ******************************************[******************************
if p.options.bmibnb.verbose>0
	disp('* -More root-node bound-propagation');   
end
p = presolveloop(p,upper);

% *************************************************************************
% Try a little more root node tigthening
% *************************************************************************
close = find(abs(p.lb - p.ub) < 1e-12);
p.lb(close) = (p.lb(close)+p.ub(close))/2;
p.ub(close) = p.lb(close);
if p.options.bmibnb.verbose>0 && p.options.bmibnb.roottight
    disp('* -Performing LP-based bound-propagation ');
end
[p,timing] = root_node_tighten(p,upper,timing);
if p.options.bmibnb.verbose>0
	disp('* -And some more root-node bound-propagation');   
end
p = propagate_bounds_from_monomials(p);
p = propagate_bounds_from_evaluations(p);
p = propagate_bounds_from_equalities(p);
p = propagate_bounds_from_monomials(p);
output = createOutputStructure(0);

% Detect complementary constraints
p.complementary = [];
complementary = find(p.lb==0 & p.ub==0 & p.variabletype'==1);
if ~isempty(complementary)
    p.complementaryvar = complementary;
    for i = 1:length(complementary)
        index = find(p.bilinears(:,1) == complementary(i));
        x = p.bilinears(index,2);
        y = p.bilinears(index,3);
        if ~((p.lb(x)==0 & p.ub(x)==0) | (p.lb(y)==0 & p.ub(y)==0))
        p.complementary = [p.complementary ;x y];
        end
    end
end

% Some code to extract bounds in complementary conditions, when the
% complementary variables hasn't been defined
for i = 1:size(p.complementary,1)
    if (isinf(p.lb(p.complementary(i,1))) |  isinf(p.ub(p.complementary(i,1))))
        if ismember(p.complementary(i,1),p.branch_variables)
            %x*y == 0 where x is unbounded. Hence x is either 0, or bounded by
            %the value obtained when y is 0
            pp=p;
            pp.ub(p.complementary(i,2))=0;
            pp = presolveloop(pp,upper);
            if pp.feasible
                p.ub(p.complementary(i,1)) = min(pp.ub(p.complementary(i,1)),p.ub(p.complementary(i,1)));
                p.lb(p.complementary(i,1)) = max(pp.lb(p.complementary(i,1)),p.lb(p.complementary(i,1)));
            else
                p.ub(p.complementary(i,1))=0;
                p.lb(p.complementary(i,1))=0;
            end
        end
    end
    if (isinf(p.lb(p.complementary(i,2))) |  isinf(p.ub(p.complementary(i,2))))
        if ismember(p.complementary(i,2),p.branch_variables)
            %x*y == 0 where y is unbounded. Hence y is either 0, or bounded by
            %the bounds obtained when x is 0
            pp=p;
            pp.ub(p.complementary(i,1))=0;
            pp = presolveloop(pp,upper);
            if pp.feasible
                p.ub(p.complementary(i,2)) = min(pp.ub(p.complementary(i,2)),p.ub(p.complementary(i,2)));
                p.lb(p.complementary(i,2)) = max(pp.lb(p.complementary(i,2)),p.lb(p.complementary(i,2)));
            else
                % This case is infeasible, hence the other term has to be zero
                p.ub(p.complementary(i,2))=0;
                p.lb(p.complementary(i,2))=0;
            end
        end
    end
end

lb = p.lb;
ub = p.ub;
if p.feasible
        
    % *******************************
    % Save time & memory
    % *******************************
    p.options.savesolverinput  = 0;
    p.options.savesolveroutput = 0;
    p.options.saveduals = 0;
    p.options.dimacs = 0;    
    
    % *******************************
    % RUN BILINEAR BRANCH & BOUND
    % *******************************
    if upper < p.options.bmibnb.target
        % Solution already good enough, no reason to start branching at all
        if p.options.bmibnb.verbose>0
            disp('* -Terminating in root as upper bound already meets target.');   
        end
        solved_nodes = 0;
        lower = -inf;
        lower_hist = -inf;
        upper_hist = upper;
        problem = 0;
        counter = p.counter;
    elseif ~isinf(upper) && isequal(p.options.bmibnb.lowertarget,-inf)
         % We only want a solution
        if p.options.bmibnb.verbose>0
            disp('* -Terminating in root as solution found and lower target is -inf');
        end
        solved_nodes = 0;
        lower = -inf;
        lower_hist = -inf;
        upper_hist = upper;
        problem = 0;
        counter = p.counter;
    else
        [x_min,solved_nodes,lower,upper,lower_hist,upper_hist,solution_hist,timing,counter,problem] = bmibnb_branch_and_bound(p,x_min,upper,timing,solution_hist);
        
        % ********************************
        % ADJUST DIAGNOSTICS
        % ********************************
        if isinf(upper) && problem == 0
            problem = 1;
        end
        if isinf(lower) & (lower<0) && problem == 0
            problem = 2;
        end
        if isinf(lower) & (lower>0) && problem == 0
            problem = 1;
        end
    end
else
    if ~isinf(upper)
        if p.options.bmibnb.verbose>0
            disp('* -Solution proven optimal already in root-node.');   
        end
        counter = p.counter;
        problem = 0;        
        solved_nodes = 0;
        lower = nan;
        lower_hist = [-inf];
        upper_hist = [upper];        
    else
        counter = p.counter;
        problem = 1;
        x_min = repmat(nan,length(p.c),1);
        solved_nodes = 0;
        lower = inf;
        lower_hist = [];
        upper_hist = [];
        solution_hist = [];
    end
end

timing.total = toc(timing.total);
if p.options.bmibnb.verbose
    disp(['* Timing: ' num2str(ceil(100*timing.uppersolve/timing.total)) '% spent in upper solver (' num2str(counter.uppersolved) ' problems solved)']);
    disp(['*         ' num2str(ceil(100*timing.lowersolve/timing.total)) '% spent in lower solver (' num2str(counter.lowersolved) ' problems solved)']);
    disp(['*         ' num2str(ceil(100*timing.lpsolve/timing.total)) '% spent in LP-based domain reduction (' num2str(counter.lpsolved) ' problems solved)']);
    disp(['*         ' num2str(ceil(100*timing.heuristics/timing.total)) '% spent in upper heuristics (' num2str(counter.heuristics) ' candidates tried)']);
end

x_min = dediagonalize(p,x_min);

output = createOutputStructure(problem);
output.solved_nodes = solved_nodes;
output.Primal        = zeros(length(p.kept),1);
output.Primal(p.kept(1:n_in))= x_min(1:n_in);
output.infostr      = yalmiperror(output.problem,'BMIBNB');
output.solvertime   = timing.total;
output.timing = timing;
output.lower = lower;
output.solveroutput.nodes = length(lower_hist);
output.solveroutput.counter = counter;
output.solveroutput.timing = timing;
output.solveroutput.lower = lower;
output.solveroutput.lower_hist = lower_hist;
output.solveroutput.upper_hist = upper_hist;
output.solveroutput.solution_hist = solution_hist;
output.extra.propagatedlb = lb;
output.extra.propagatedub = ub;


function p = presolveloop(p,upper)
i = 0;
goon = 1;
while goon && any(abs(p.ub(p.branch_variables)-p.lb(p.branch_variables))>p.options.bmibnb.vartol)
    start = [p.lb;p.ub];    
    p = propagate_bounds_from_upper(p,upper);    
    p = propagate_bounds_from_complementary(p);
    p = propagate_bounds_from_monomials(p);  
    p = propagate_bounds_from_equalities(p);       
    p = propagate_bounds_from_monomials(p);
    p = propagate_bounds_from_monomials(p);        
    p = propagate_bounds_from_evaluations(p);    
    goon = (norm(start-[p.lb;p.ub],'inf') > 1e-2) & i < 8;    
    i = i+1;
end


function x = dediagonalize(p,x);
if isempty(p.V)
    return
end
y = p.V*x(1:length(x)/2);
x = zeros(p.original_n,1);
x(p.original_linear) = y;