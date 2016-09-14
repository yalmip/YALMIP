function output = bmibnb(p)
%BMIBNB Global solver based on branch-and-bound
%
% BMIBNB is never called by the user directly, but is called by YALMIP from
% SOLVESDP, by choosing the solver tag 'bmibnb' in sdpsettings
%
% The behaviour of BMIBNB can be altered using the fields in the field
% 'bmibnb' in SDPSETTINGS
%
% bmibnb.lowersolver    - Solver for lower bound [solver tag ('')]
% bmibnb.uppersolver    - Solver for upper bound [solver tag ('')]
% bmibnb.lpsolver       - Solver for LP bound tightening [solver tag ('')]
% bmibnb.branchmethod   - Branch strategy ['maxvol' | 'best' ('best')]
% bmibnb.branchrule     - Branch position ['omega' | 'bisect' ('omega')]
% bmibnb.lpreduce       - Improve variable bounds using LP [ real [0,1] (0 means no reduction, 1 means all variables)
% bmibnb.lowrank        - Partition variables into two disjoint sets and branch on smallest [ 0|1 (0)]
% bmibnb.target         - Exit if upper bound<target [double (-inf)]
% bmibnb.roottight      - Improve variable bounds in root using full problem [ 0|1 (1)]
% bmibnb.vartol         - Cut tree when x_U-x_L < vartol on all branching variables
% bmibnb.relgaptol      - Tolerance on relative objective error (UPPER-LOWER)/(1+|UPPER|) [real (0.01)]
% bmibnb.absgaptol      - Tolerance on objective error (UPPER-LOWER) [real (0.01)]
% bmibnb.pdtol          - A number is declared non-negative if larger than...[ double (-1e-6)]
% bmibnb.eqtol          - A number is declared zero if abs(x) smaller than...[ double (1e-6)]
% bmibnb.maxiter        - Maximum number nodes [int (100)]
% bmibnb.maxtime        - Maximum CPU time (sec.) [int (3600)]

% Author Johan Löfberg

bnbsolvertime = clock;
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

% CPLEX handles options insanely slow, so we remove all default values in
% case cplex will be used
p.options = prunecplexoptions(p.options);

timing.total = tic;
timing.uppersolve = 0;
timing.lowersolve = 0;
timing.domainreduce = 0;
if ~isempty(p.F_struc)
    if any(isnan(p.F_struc) | isinf(p.F_struc))
        output = yalmip_default_output;
        output.problem = 1;
        output.solved_nodes = 0;
        output.Primal       = zeros(length(p.c),1);
        output.infostr      = yalmiperror(output.problem,'BNB');
        output.solvertime   = 0;
        output.lower = inf;
        return
    end
end

% *************************************************************************
% Assume feasible (this property can be changed in the presolve codes
% *************************************************************************
p.feasible = 1;

% *************************************************************************
% Ensure that bond propagation is performed
% *************************************************************************
p.changedbounds = 1;

% *************************************************************************
% Initialize some diagnostic counters
% *************************************************************************
p.counter.lpsolved = 0;
p.counter.lowersolved = 0;
p.counter.uppersolved = 0;
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
p = presolve_bounds_from_domains(p);
p = presolve_bounds_from_modelbounds(p);
%p = presolve_eliminatelinearratios(p);
p = presolve_bounds_from_quadratics(p);
p = update_eval_bounds(p);
p = update_monomial_bounds(p);


% *************************************************************************
% Sort equalities in order to improve future bound propagation
% *************************************************************************
p = presolve_sortrows(p);

% *************************************************************************
% Improve the bounds by performing some standard propagation
% *************************************************************************
p = propagate_bounds_from_equalities(p); 
p = update_eval_bounds(p);
p = update_sumsepquad_bounds(p);
p = update_monomial_bounds(p);
p = presolve_bounds_from_inequalities(p);
p = update_eval_bounds(p);
p = update_monomial_bounds(p);
p = presolve_quadratic_psdbound(p);
% *************************************************************************
% For quadratic nonconvex programming over linear constraints, we
% diagonalize the problem to obtain less number of bilinear terms. Not
% theoretically any reason to do so, but practical performance is better
% *************************************************************************
p = diagonalize_quadratic_program(p);
if p.diagonalized
    p = compile_nonlinear_table(p);
end

% *************************************************************************
% Try to generate a feasible solution, by using avialable x0 (if usex0=1),
% or by trying the zero solution, or my trying the (lb+ub)/2 solution. Note
% that we do this here b4efore the problem possibly is bilinearized, thus
% avoiding to introduce possibly complicating bilinear constraints
% *************************************************************************
[p,x_min,upper] = initializesolution(p);

% *************************************************************************
% If the upper bound solver is capable of solving the original problem,
% without bilinearizing it first, we might try to get a local solution,
% possibly based on our initial crude solution
% *************************************************************************
if solver_can_solve(p.solver.uppersolver,p) & any(p.variabletype>2)
    p.high_monom_model = [];
    p = build_recursive_scheme(p);
    p = compile_nonlinear_table(p);
    p = preprocess_bilinear_bounds(p);
    p = update_eval_bounds(p);
    p = updateboundsfromupper(p,upper);
    p = propagate_bounds_from_complementary(p);
    p = updatemonomialbounds(p);
    p = propagate_bounds_from_equalities(p);
    p = updatemonomialbounds(p);   
    p = updatemonomialbounds(p);    
    p = update_eval_bounds(p);    
    [upper,p.x0,info_text,numglobals,timing] = solve_upper_in_node(p,p,x_min,upper,x_min,p.solver.uppersolver.call,'',0,timing);
    p.counter.uppersolved = p.counter.uppersolved + 1;
    if numglobals > 0
        x_min = p.x0;
    end
end
if isempty(p.x0)
    p.x0 = zeros(length(p.c),1);
end

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
p = propagatequadratics(p);
    

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

% *************************************************************************
% Tighten bounds (might be useful after bilinearization?)
% *************************************************************************
p = preprocess_bilinear_bounds(p);
p = update_eval_bounds(p);

% *************************************************************************
% Now reduce the branch variables by removing bilinear terms that only have
% been introduced due to a convex quadratic objective
% *************************************************************************
p = reduce_bilinear_branching_variables(p);

% *************************************************************************
% Simple pre-solve loop. The loop is needed for variables defined as w =
% x*y, x = t*u,y=..db
% ******************************************[******************************
p = presolveloop(p,upper);

% *************************************************************************
% Try a little more root node tigthening
% *************************************************************************
close = find(abs(p.lb - p.ub) < 1e-12);
p.lb(close) = (p.lb(close)+p.ub(close))/2;
p.ub(close) = p.lb(close);
p = root_node_tighten(p,upper);
p = updatemonomialbounds(p);
p = update_eval_bounds(p);
p = propagate_bounds_from_equalities(p);
p = updatemonomialbounds(p);
output = yalmip_default_output;

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
    % Bounded domain?
    % *******************************
    if ~isempty(p.branch_variables)
        if any(isinf(p.lb(p.branch_variables))) | any(isinf(p.ub(p.branch_variables)))
            bad = find(isinf(p.lb(p.branch_variables)) | isinf(p.ub(p.branch_variables)));
            if ~all(ismember(bad,p.complementary))
                output = yalmip_default_output;
                output.Primal  = x_min;
                output.problem = -6;
                output.infostr = yalmiperror(-6);
                output.solved_nodes = 0;
                return
            end
        end
    end

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
    [x_min,solved_nodes,lower,upper,lower_hist,upper_hist,timing,counter] = branch_and_bound(p,x_min,upper,timing);
       
    % ********************************
    % CREATE SOLUTION AND DIAGNOSTICS
    % ********************************
    problem = 0;
    if isinf(upper)
        problem = 1;
    end
    if isinf(lower) & (lower<0)
        problem = 2;
    end
    if isinf(lower) & (lower>0)
        problem = 1;
    end    
    if solved_nodes == p.options.bmibnb.maxiter
        problem = 3;
    end
else
    counter = p.counter;
    problem = 1;
    x_min = repmat(nan,length(p.c),1);
    solved_nodes = 0;
    lower = inf;
    lower_hist = [];
    upper_hist = [];
end

timing.total = toc(timing.total);
if p.options.bmibnb.verbose
    disp(['* Timing: ' num2str(ceil(100*timing.uppersolve/timing.total)) '% spent in upper solver (' num2str(counter.uppersolved) ' problems solved)']);
    disp(['*         ' num2str(ceil(100*timing.lowersolve/timing.total)) '% spent in lower solver (' num2str(counter.lowersolved) ' problems solved)']);
    disp(['*         ' num2str(ceil(100*timing.domainreduce/timing.total)) '% spent in LP-based domain reduction (' num2str(counter.lpsolved) ' problems solved)']);
end

x_min = dediagonalize(p,x_min);

output = yalmip_default_output;
output.problem = problem;
output.solved_nodes = solved_nodes;
output.Primal        = zeros(length(p.kept),1);
output.Primal(p.kept(1:n_in))= x_min(1:n_in);
output.infostr      = yalmiperror(output.problem,'BMIBNB');
output.solvertime   = etime(clock,bnbsolvertime);
output.timing = timing;
output.lower = lower;
output.solveroutput.lower = lower;
output.solveroutput.lower_hist = lower_hist;
output.solveroutput.upper_hist = upper_hist;
output.extra.propagatedlb = lb;
output.extra.propagatedub = ub;

function pnew = diagonalize_quadratic_program(p);
pnew = p;
pnew.V = [];
pnew.diagonalized = 0;

% No quadratic terms
if all(p.variabletype == 0)
    return
end

% Any polynomial terms or simple linear by some reason
if any(p.variabletype > 2) & ~all(p.variabletype == 0) | p.options.bmibnb.diagonalize==0
    return
end

if ~isempty(p.evalVariables)
    return
end

if ~isempty(p.binary_variables) | ~isempty(p.integer_variables)
    return
end


nonlinear = find(p.variabletype > 0);
linear = find(p.variabletype == 0);
if ~isempty(p.F_struc)
    % Nonlinear terms in constraints
    if nnz(p.F_struc(:,1 + nonlinear))> 0
        return
    end
end

if ~isempty(p.lb)
    if ~all(isinf(p.lb(nonlinear)))
        return
    end
end
if ~isempty(p.ub)
    if ~all(isinf(p.ub(nonlinear)))
        return
    end
end

% Find quadratic and linear terms
used_in_c = find(p.c);
quadraticterms = used_in_c(find(ismember(used_in_c,nonlinear)));
Q = zeros(length(p.c),length(p.c));
if ~isempty(quadraticterms)
    usedinquadratic = zeros(1,length(p.c));
    for i = 1:length(quadraticterms)
        Qij = p.c(quadraticterms(i));
        power_index = find(p.monomtable(quadraticterms(i),:));
        if length(power_index) == 1
            Q(power_index,power_index) = Qij;
        else
            Q(power_index(1),power_index(2)) = Qij/2;
            Q(power_index(2),power_index(1)) = Qij/2;
        end       
    end
end
Qlin = Q(linear,linear);
clin = p.c(linear);

% Decompose Q
[V,D] = eig(full(Qlin));
V = real(V);
D = real(D);
V(abs(V)<1e-11) = 0;
D(abs(D)<1e-11) = 0;
lb = p.lb(linear);
ub = p.ub(linear);
Z = V';
newub = inf(length(lb),1);
newlb = -inf(length(lb),1);
for i = 1:length(lb)
    z = Z(i,:);
    neg = find(z<0);
    pos = find(z>0);
    if ~isempty(neg)
        newub(i,1) = z(neg)*lb(neg);
        newlb(i,1) = z(neg)*ub(neg);
    end
%    if ~isempty(pos)
        newub(i,1) = z(pos)*ub(pos);
        newlb(i,1) = z(pos)*lb(pos);        
%    end
%  newub(i,1) = z(z>0)*ub(z>0) + z(z<0)*lb(z<0);
%        newlb(i,1) = z(z>0)*lb(z>0) + z(z<0)*ub(z<0);  
end


% Create new problem
clin = V'*clin;
%A = A*V;

n = length(linear);
pnew.original_linear = linear;
pnew.original_n = length(p.c);
pnew.V = V;
pnew.c = [clin;diag(D)];
pnew.Q = spalloc(2*n,2*n,0);

% find constraint polytope
if size(p.F_struc,1)>0
    A = -p.F_struc(:,1 + linear);
    b = p.F_struc(:,1);
    pnew.F_struc = [b -A*V zeros(length(b),n)];
end


pnew.variabletype = [zeros(1,n) ones(1,n)*2];
pnew.monomtable = [eye(n);2*eye(n)];
pnew.monomtable(2*n,2*n) = 0;
pnew.lb =-inf(2*n,1);
pnew.ub = inf(2*n,1);
pnew.lb(1:n) = newlb;
pnew.ub(1:n) = newub;
if length(pnew.x0)>0
    pnew.x0 = V*p.x0;
end
pnew.diagonalized = 1;

function x = dediagonalize(p,x);
if isempty(p.V)
    return
end
y = p.V*x(1:length(x)/2);
x = zeros(p.original_n,1);
x(p.original_linear) = y;


function p = presolveloop(p,upper)
i = 0;
goon = 1;
while goon
    start = [p.lb;p.ub];
    i = i+1;
    % 2
    p = updateboundsfromupper(p,upper);
    % 3
    p = propagate_bounds_from_complementary(p);
    p = updatemonomialbounds(p);  
    p = propagate_bounds_from_equalities(p);       
    p = updatemonomialbounds(p);
    p = updatemonomialbounds(p);        
    p = update_eval_bounds(p);
    i = i+1;
    goon = (norm(start-[p.lb;p.ub],'inf') > 1e-2) & i < 15;
    start = [p.lb;p.ub];
end


