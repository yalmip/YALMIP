function output = bnb(p)
%BNB          General branch-and-bound scheme for (primarily) conic programs
%
% BNB applies a branch-and-bound scheme to solve mixed-integer
% convex programs, in particular linear mixed-integer semidefinite
% programs
%
% BNB is never called by the user directly, but is called by
% YALMIP from SOLVESDP, by choosing the solver tag 'bnb' in sdpsettings.
%
% BNB is used if no other mixed-integer solver is found, and
% is only meant to be used for mixed-integer SDP, or maybe general
% nonlinear mixed-integer problems (as there are much better solvers
% available for standard LP/QP/SOCP models) 
%
% The behaviour of BNB can be altered using the fields in the field 'bnb'
% in SDPSETTINGS (although defaults are recommended)
%
% solver           Solver for the relaxed problems (standard solver tag, see SDPSETTINGS)
%
% maxiter          Maximum number of nodes explored
%
% maxtime          Maximum time allowed
%
% inttol           Tolerance for declaring a variable as integer
%
% feastol          Tolerance for declaring constraints as feasible
%
% gaptol           Exit when (upper bound-lower bound)/(1e-3+abs(lower bound)) < gaptol
%
% round            Round variables smaller than bnb.inttol
%
% plot             Plot the upper and lower bound, stack size, and
%                  histogram of lower bounds in stack
%
% bnb.branchrule   Deceides on what variable to branch
%                   'max'     : Variable furthest away from being integer
%                   'min'     : Variable closest to be being integer
%                   'first'   : First variable (lowest variable index in YALMIP)
%                   'last'    : Last variable (highest variable index in YALMIP)
%                   'weight'  : See manual
%
% bnb.method       Branching strategy
%                   'depth'   : Depth first
%                   'breadth' : Breadth first
%                   'best'    : Expand branch with lowest lower bound
%                   'depthX'  : Depth until integer solution found, then X (e.g 'depthbest')
%
%
% See also OPTIMIZE, BINVAR, INTVAR, BINARY, INTEGER

% ********************************
%% INITIALIZE DIAGNOSTICS IN YALMIP
% ********************************
bnbsolvertime = clock;
showprogress('Branch and bound started',p.options.showprogress);

% ********************************
%% Remove options if none has been changed
%%  Improves peroformance when calling solver many times
% ********************************
p.options = pruneOptions(p.options);

% ********************************
%% We might have a GP : pre-calc
% ********************************
p.nonlinear = find(~(sum(p.monomtable~=0,2)==1 & sum(p.monomtable,2)==1));
p.nonlinear = union(p.nonlinear,p.evalVariables);

% ********************************
% This field is only used in bmibnb, which uses the same sub-functions as
% bnb
% ********************************
p.high_monom_model = [];

% ********************************
%% Define infinite bounds
% ********************************
if isempty(p.ub)
    p.ub = repmat(inf,length(p.c),1);
end
if isempty(p.lb)
    p.lb = repmat(-inf,length(p.c),1);
end

% ********************************
%% Extract bounds from model
% ********************************
p = extractBounds(p);

% ********************************
%% ADD CONSTRAINTS 0<x<1 FOR BINARY
% ********************************
if ~isempty(p.binary_variables)
    p.ub(p.binary_variables) =  min(p.ub(p.binary_variables),1);
    p.lb(p.binary_variables) =  max(p.lb(p.binary_variables),0);     
end

p = update_integer_bounds(p);

if ~isempty(p.semicont_variables)
    redundant = find(p.lb<=0 & p.ub>=0);
    p.semicont_variables = setdiff(p.semicont_variables,redundant);
    % Now relax the model and generate hull including 0
    p.semibounds.lb = p.lb(p.semicont_variables);
    p.semibounds.ub = p.ub(p.semicont_variables);
    p.lb(p.semicont_variables) = min(p.lb(p.semicont_variables),0);
    p.ub(p.semicont_variables) = max(p.ub(p.semicont_variables),0);
end

% Could be some nonlinear terms (although these problems are recommended to
% be solved using BMIBNB
p = compile_nonlinear_table(p);
p = updatemonomialbounds(p);

% % *******************************
% %% PRE-SOLVE (nothing fancy coded)
% % *******************************
p = simplePresolve(p);
p = propagate_bounds_from_equalities(p);
pss = [];
if isempty(p.nonlinear)
    if p.K.f>0
        Aeq = -p.F_struc(1:p.K.f,2:end);
        beq = p.F_struc(1:p.K.f,1);
        A = [Aeq;-Aeq];
        b = [beq;-beq];
        [p.lb,p.ub,redundant,pss] = tightenbounds(A,b,p.lb,p.ub,p.integer_variables,p.binary_variables,ones(length(p.lb),1));
    end
    pss=[];
    if p.K.l>0
        A = -p.F_struc(1+p.K.f:p.K.f+p.K.l,2:end);
        b = p.F_struc(1+p.K.f:p.K.f+p.K.l,1);
        [p.lb,p.ub,redundant,pss] = tightenbounds(A,b,p.lb,p.ub,p.integer_variables,p.binary_variables,ones(length(p.lb),1));
        if length(redundant)>0
            pss.AL0A(redundant,:)=[];
            pss.AG0A(redundant,:)=[];
            p.F_struc(p.K.f+redundant,:)=[];
            p.K.l = p.K.l - length(redundant);
        end
    end
end

% Silly redundancy
p = updatemonomialbounds(p);
p = propagate_bounds_from_equalities(p);
if p.K.l > 0
    b = p.F_struc(1+p.K.f:p.K.l+p.K.f,1);
    A = -p.F_struc(1+p.K.f:p.K.l+p.K.f,2:end);
    redundant = find(((A>0).*A*(p.ub-p.lb) - (b-A*p.lb) <= 0));
    if ~isempty(redundant)
        p.F_struc(p.K.f + redundant,:) = [];
        p.K.l = p.K.l - length(redundant);
    end
end

% *******************************
%% Display logics
% 0 : Silent
% 1 : Display branching
% 2 : Display node solver prints
% *******************************
if p.options.verbose ~= fix(p.options.verbose)
    p.options.print_interval = ceil(1/p.options.verbose);
    p.options.verbose = ceil(p.options.verbose);
else
    p.options.print_interval = 1;
end
switch max(min(p.options.verbose,3),0)
    case 0
        p.options.bnb.verbose = 0;
    case 1
        p.options.bnb.verbose = 1;
        p.options.verbose = 0;
    case 2
        p.options.bnb.verbose = 2;
        p.options.verbose = 0;
    case 3
        p.options.bnb.verbose = 2;
        p.options.verbose = 1;
    otherwise
        p.options.bnb.verbose = 0;
        p.options.verbose = 0;
end

% *******************************
%% Figure out the weights if any
% *******************************
try % Probably buggy first version...
    if ~isempty(p.options.bnb.weight)
        weightvar = p.options.bnb.weight;
        if isa(weightvar,'sdpvar')
            if (prod(size(weightvar)) == 1)
                weight = ones(length(p.c),1);
                for i = 1:length(p.c)
                    weight(i,1) = full(getbasematrix(weightvar,p.used_variables(i)));
                end
                p.weight = weight;
            else
                error('Weight should be an SDPVAR scalar');
            end
        else
            error('Weight should be an SDPVAR scalar');
        end
    else
        p.weight = ones(length(p.c),1);
    end
catch
    disp('Something wrong with weights. Please report bug');
    p.weight = ones(length(p.c),1);
end

% *******************************
%% START BRANCHING
% *******************************
setuptime = etime(clock,bnbsolvertime);
bnbsolvertime = clock;
[x_min,solved_nodes,lower,upper,profile,diagnostics] = branch_and_bound(p,pss);
bnbsolvertime =  etime(clock,bnbsolvertime);
output.solvertime = setuptime + bnbsolvertime;

% **********************************
%% CREATE SOLUTION
% **********************************
if diagnostics == -4
    output.problem = -4;
else
    output.problem = 0;
    if isinf(upper)
        output.problem = 1;
    end
    if isinf(-lower)
        output.problem = 2;
    end
    if solved_nodes == p.options.bnb.maxiter
        output.problem = 3;
    end
    if bnbsolvertime > p.options.bnb.maxtime
        output.problem = 3;
    end
end
output.solved_nodes = solved_nodes;
output.Primal      = x_min;
output.Dual        = [];
output.Slack       = [];
if  output.problem == -4
    output.infostr      = yalmiperror(output.problem,[p.solver.lower.tag '-' p.solver.lower.version]);
else
    output.infostr      = yalmiperror(output.problem,'BNB');
end
output.solverinput  = 0;
if p.options.savesolveroutput
    output.solveroutput.setuptime = setuptime;
    output.solveroutput.localsolvertime = profile.local_solver_time;
    output.solveroutput.branchingtime = bnbsolvertime;
    output.solveroutput.solved_nodes = solved_nodes;
    output.solveroutput.lower = lower;
    output.solveroutput.upper = upper;
else
    output.solveroutput =[];
end
%% --

function [x_min,solved_nodes,lower,upper,profile,diagnostics] = branch_and_bound(p,pss)

% *******************************
% We don't need this
% *******************************
p.options.savesolveroutput = 0;
p.options.saveduals = 0;
p.options.dimacs = 0;
diagnostics = 0;
bnbsolvertime = clock;

% *******************************
% Tracking performance etc
% *******************************
profile.local_solver_time = 0;

% *************************************************************************
% We save this to re-use some stuff in fmincon
% *************************************************************************
p.options.savesolverinput = 1;

% *******************************
%% SET-UP ROOT PROBLEM
% *******************************
p.depth = 0;
p.lower = NaN;
% Does the user want to create his own initial guess
if p.options.usex0
    [x_min,upper] = initializesolution(p);
    if isinf(upper)
        % Try to initialize to lowerbound+upperbound. fmincon really
        % doesn't like zero initial guess, despite having bounds available
        x_min = zeros(length(p.c),1);
        violates_finite_bounds = ((x_min < p.lb) | (x_min < p.ub));
        violates_finite_bounds = find(violates_finite_bounds & ~isinf(p.lb) & ~isinf(p.ub));
        x_min(violates_finite_bounds) = (p.lb(violates_finite_bounds) + p.ub(violates_finite_bounds))/2;
        x_min = setnonlinearvariables(p,x_min);
    end
    p.x0 = x_min;
else
    upper = inf;
    x_min = zeros(length(p.c),1);
    violates_finite_bounds = ((x_min < p.lb) | (x_min < p.ub));
    violates_finite_bounds = find(violates_finite_bounds & ~isinf(p.lb) & ~isinf(p.ub));
    x_min(violates_finite_bounds) = (p.lb(violates_finite_bounds) + p.ub(violates_finite_bounds))/2;
    x_min = setnonlinearvariables(p,x_min);
    p.x0 = x_min;
end


% *******************************
%% Global stuff
% *******************************
lower = NaN;
stack = stackCreate;

% *******************************
%% Create function handle to solver
% *******************************
lowersolver = p.solver.lower.call;
uppersolver = p.options.bnb.uppersolver;

p.corig = p.c;

% *******************************
%% INVARIANT PROBLEM DATA
% *******************************
c = p.corig;
Q = p.Q;
f = p.f;
integer_variables = p.integer_variables;
solved_nodes = 0;

semicont_variables = p.semicont_variables;

gap = inf;
node = 1;

if p.options.bnb.presolve
    savec = p.c;
    saveQ = p.Q;
    p.Q = p.Q*0;
    
    n = length(p.c);
    saveBinary = p.binary_variables;
    saveInteger = p.integer_variables;
    p.binary_variables = [];
    p.integer_variables = [];;
    
    for i = 1:length(c)
        p.c = eyev(n,i);
        output = feval(lowersolver,p);
        if output.problem == 0
            p.lb(i) = max(p.lb(i),output.Primal(i)-1e-3);
        end
        p.c = -eyev(n,i);
        output = feval(lowersolver,p);
        if output.problem == 0
            p.ub(i) = min(p.ub(i),output.Primal(i)+1e-3);
        end
        p.lb(saveBinary) = ceil(p.lb(saveBinary)-1e-3);
        p.ub(saveBinary) = floor(p.ub(saveBinary)+1e-3);
    end
    p.binary_variables = saveBinary;
    p.integer_variables = saveInteger;
    
    p.Q = saveQ;
    p.c = savec;
end

% ************************************************
% Some hacks to speed up solver calls
% Only track solver-time if user wants profile
% ************************************************
p.getsolvertime = p.options.bnb.profile;

% *******************************
%% DISPLAY HEADER
% *******************************
originalDiscrete = [p.integer_variables(:);p.binary_variables(:)];
originalBinary   = p.binary_variables(:);

if nnz(Q)==0 & (nnz(p.c-fix(p.c))==0) & isequal(p.K.m,0)
    can_use_ceil_lower = all(ismember(find(p.c),originalDiscrete));
else
    can_use_ceil_lower = 0;
end

if p.options.bnb.verbose
    
    pc = p.problemclass;
    non_convex_obj = pc.objective.quadratic.nonconvex | pc.objective.polynomial;
    non_convex_constraint =  pc.constraint.equalities.quadratic | pc.constraint.inequalities.elementwise.quadratic.nonconvex;
    non_convex_constraint =  non_convex_constraint | pc.constraint.equalities.polynomial | pc.constraint.inequalities.elementwise.polynomial;
    
    possiblynonconvex = non_convex_obj | non_convex_constraint;
    if ~isequal(p.solver.lower.version,'')
        p.solver.lower.tag = [p.solver.lower.tag '-' p.solver.lower.version];
    end
    
    disp('* Starting YALMIP integer branch & bound.');
    disp(['* Lower solver   : ' p.solver.lower.tag]);
    disp(['* Upper solver   : ' p.options.bnb.uppersolver]);
    disp(['* Max time       : ' num2str(p.options.bnb.maxtime)]);
    disp(['* Max iterations : ' num2str(p.options.bnb.maxiter)]);
    
    if possiblynonconvex & p.options.warning
        disp(' ');
        disp('Warning : The continuous relaxation may be nonconvex. This means ');
        disp('that the branching process is not guaranteed to find a');
        disp('globally optimal solution, since the lower bound can be');
        disp('invalid. Hence, do not trust the bound or the gap...')
        
    end
end
if p.options.bnb.verbose;            disp(' Node       Upper       Gap(%)     Lower     Open   Elapsed time');end;

if nnz(Q)==0 & nnz(c)==1 & isequal(p.K.m,0)
    p.simplecost = 1;
else
    p.simplecost = 0;
end

poriginal = p;
p.cuts = [];

p = detectSOS(p);
p = detectAtMost(p);
poriginal.atmost = p.atmost;

pid = 0;
lowerhist = [];
upperhist = [];
stacksizehist = [];
p.fixedvariable = [];
p.fixdir = '';
lastUpper = upper;

oldp = p;

if length(p.integer_variables) == length(p.c)
    p.all_integers = 1;
else
    p.all_integers = 0;
end
p.noninteger_variables = setdiff(1:length(p.c),[p.integer_variables p.binary_variables p.semicont_variables]);
poriginal.noninteger_variables = p.noninteger_variables;

p = addImpliedSDP(p);

% Resuse some code from cutsdp to add simple cuts required for SDP
% feasibility for problems with some trivial symmetries
% TODO: Clean up, refactor, generalize
if p.K.f == 0 % Still to lazy to fix last insertion
    top = 1 + p.K.l + p.K.f + sum(p.K.q);
    p.F_struc = p.F_struc';
    for i = 1:length(p.K.s)
        p.semidefinite{i}.F_struc = p.F_struc(:,top:top+p.K.s(i)^2-1)';
        top = top + p.K.s(i)^2;
    end
    p.F_struc = p.F_struc';
    p.sdpsymmetry = [];
    p = detect3x3SymmetryGroups(p);
    pp = p;pp.F_struc = [];pp.K.l = 0;pp.K.f = 0;pp.K.s = 0;
    [p,pp] = addSymmetryCuts(p,pp);
    p = addLinearCut(p,pp.F_struc);
    p.semidefinite=[];
end

feasibilityHistory = [];
% Not used, deleted feature...
aggresiveprune = 0;
% Save of all optimal solutions
allSolutions = [];
sosgroups = [];
sosvariables = [];
while ~isempty(node) & (etime(clock,bnbsolvertime) < p.options.bnb.maxtime) & (solved_nodes < p.options.bnb.maxiter) & (isinf(lower) | gap>p.options.bnb.gaptol)
    
    % ********************************************
    % BINARY VARIABLES ARE FIXED ALONG THE PROCESS
    % ********************************************
    binary_variables = p.binary_variables;
    
    % ********************************************
    % SO ARE SEMI VARIABLES
    % ********************************************
    semicont_variables  = p.semicont_variables;
    
    % ********************************************
    % ASSUME THAT WE WON'T FATHOME
    % ********************************************
    keep_digging = 1;
    message = '';
    
    % *************************************
    % SOLVE NODE PROBLEM
    % *************************************
    if any(p.ub<p.lb - 1e-12)
        x = zeros(length(p.c),1);
        output.Primal = x;
        output.problem=1;
    else
        p.x_min = x_min;
        relaxed_p = p;
        relaxed_p.integer_variables = [];
        relaxed_p.binary_variables = [];
        relaxed_p.semicont_variables = [];
        relaxed_p.ub(p.ub<p.lb) = relaxed_p.lb(p.ub<p.lb);
        
        % Solve node relaxation     
        output = bnb_solvelower(lowersolver,relaxed_p,upper,lower,x_min,aggresiveprune,allSolutions);
        if (output.problem == 12 || output.problem == 2) && ~isinf(p.lower)
            output.problem = 1;
        end
       
        if p.options.bnb.profile
            profile.local_solver_time  = profile.local_solver_time + output.solvertime;
        end
        
        % A bit crappy code to exploit computations that were done in the
        % call to fmincon...
        if isfield(output,'solverinput')
            if isfield(output.solverinput,'model')
                if isfield(output.solverinput.model,'fastdiff')
                    p.fastdiff = output.solverinput.model.fastdiff;
                end
            end
        end
                
        if output.problem == -4
            diagnostics = -4;
            x = nan+zeros(length(p.lb),1);
        else
            if isempty(output.Primal)
                output.Primal = zeros(length(p.c),1);
            end
            x  = setnonlinearvariables(p,output.Primal);
            if(p.K.l>0) & any(p.F_struc(p.K.f+1:p.K.f+p.K.l,:)*[1;x]<-1e-5)
                output.problem = 1;
            elseif output.problem == 5 & ~checkfeasiblefast(p,x,p.options.bnb.feastol)
                output.problem = 1;
            end
        end
    end
    
    solved_nodes = solved_nodes+1;
    
    % **************************************
    % THIS WILL BE INTIAL GUESS FOR CHILDREN
    % **************************************
    p.x0 = x;
    
    % *************************************
    % ANY INTEGERS? ROUND?
    % *************************************
    non_integer_binary = abs(x(binary_variables)-round(x(binary_variables)))>p.options.bnb.inttol;
    non_integer_integer = abs(x(integer_variables)-round(x(integer_variables)))>p.options.bnb.inttol;
    if p.options.bnb.round
        x(binary_variables(~non_integer_binary))   = round(x(binary_variables(~non_integer_binary)));
        x(integer_variables(~non_integer_integer)) = round(x(integer_variables(~non_integer_integer)));
    end
    non_integer_binary = find(non_integer_binary);
    non_integer_integer = find(non_integer_integer);
    if isempty(p.semicont_variables)
        non_semivar_semivar=[];
    else
        non_semivar_semivar = find(~(abs(x(p.semicont_variables))<p.options.bnb.inttol | (x(p.semicont_variables)>p.semibounds.lb & x(p.semicont_variables)<=p.semibounds.ub)));
    end
    
    x  = setnonlinearvariables(p,x);
    
    TotalIntegerInfeas = sum(abs(round(x(non_integer_integer))-x(non_integer_integer)));
    TotalBinaryInfeas = sum(abs(round(x(non_integer_binary))-x(non_integer_binary)));
    
    % *************************************
    % NODE HEURISTICS (NOTHING CODED)
    % *************************************
    should_be_tight = find([p.lb == p.ub]);
    if ~isempty(should_be_tight)
        % FIX for problems that only report numerical problems but violate
        % binary
        if max(abs(p.lb(should_be_tight)-x(should_be_tight)))>p.options.bnb.inttol
            output.problem = 1;
        end
    end
    
    if output.problem==0 | output.problem==3 | output.problem==4
        cost = computecost(f,c,Q,x,p);
        
        if output.problem~=1
            if isnan(lower)
                lower = cost;
            end
            
            if cost <= upper & ~(isempty(non_integer_binary) & isempty(non_integer_integer) & isempty(non_semivar_semivar))
                poriginal.upper = upper;
                poriginal.lower = lower;
                [upper1,x_min1] = feval(uppersolver,poriginal,output,p);               
                if upper1 < upper
                    x_min = x_min1;
                    allSolutions = x_min;
                    upper = upper1;
                    [stack,stacklower] = prune(stack,upper,p.options,solved_nodes,p,allSolutions);
                    lower = min(lower,stacklower);                   
                elseif ~isinf(upper1) && upper1 == upper && norm(x_min-x_min1) > 1e-4;
                    % Yet another solution with same value
                     allSolutions = [allSolutions x_min1];                    
                end
            elseif isempty(non_integer_binary) && isempty(non_integer_integer) && isempty(non_semivar_semivar)
            end
        end
    end
        
    p = adaptivestrategy(p,upper,solved_nodes);
    
    % *************************************
    % CHECK FATHOMING POSSIBILITIES
    % *************************************
    feasible = 1;
    
    switch output.problem
        case {-1,4}
            % Solver behaved weird. Make sure we continue digging
            keep_digging = 1;
            feasible = 1;
            cost = lower;
            x = p.lb + (p.ub-p.lb)*(1/pi);
        case 0
            if can_use_ceil_lower
                lower = ceil(lower-1e-8);
            end
        case {1,12,-4,22}
            keep_digging = 0;
            cost = inf;
            feasible = 0;
        case 2
            cost = -inf;
        otherwise
            % This part has to be much more robust because this could be a
            % numerical problem leading to a non-trustworthy solution
            cost = f+c'*x+x'*Q*x;
    end
    
    % **************************************
    % YAHOO! INTEGER SOLUTION FOUND
    % **************************************
    if isempty(non_integer_binary) & isempty(non_integer_integer)  & isempty(non_semivar_semivar) & ~(output.problem == -1) &  ~(output.problem == 4)
        if (cost<upper) & feasible
            x_min = x;
            upper = cost;
            allSolutions = x_min;
            [stack,lower] = prune(stack,upper,p.options,solved_nodes,p,allSolutions);
        end
        p = adaptivestrategy(p,upper,solved_nodes);
        keep_digging = 0;
    end
    
    % **************************************
    % Stop digging if it won't give sufficient improvement anyway
    % **************************************
    if cost>upper*(1-p.options.bnb.gaptol)
        keep_digging = 0;
    end
    
    feasibilityHistory(end+1) = feasible;
    % **********************************
    % CONTINUE SPLITTING?
    % **********************************
    if keep_digging & (cost<upper)
        
        if solved_nodes == 1
            RootNodeInfeas =  TotalIntegerInfeas+TotalBinaryInfeas;
            RootNodeCost = cost;
        end
        
        % **********************************
        % BRANCH VARIABLE
        % **********************************
        [index,whatsplit,globalindex] = branchvariable(x,integer_variables,binary_variables,p.options,x_min,[],p);
        
        % **********************************
        % CREATE NEW PROBLEMS
        % **********************************
        p0_feasible = 1;
        p1_feasible = 1;
        
        switch whatsplit
            case 'binary'
                [p0,p1,index] = binarysplit(p,x,index,cost,[],sosgroups,sosvariables);
                
            case 'integer'
                [p0,p1] = integersplit(p,x,index,cost,x_min);
                
            case 'semi'
                [p0,p1] = semisplit(p,x,index,cost,x_min);
                
            case 'sos1'
                [p0,p1] = sos1split(p,x,index,cost,x_min);
                                        
            otherwise
        end
                
        node1.lb = p1.lb;
        node1.ub = p1.ub;
        node1.depth = p1.depth;
        node1.lower = p1.lower;
        node1.fixedvariable = globalindex;
        node1.fixdir = 'up';
        node1.TotalIntegerInfeas = TotalIntegerInfeas;
        node1.TotalBinaryInfeas = TotalBinaryInfeas;
        node1.IntInfeas = 1-(x(globalindex)-floor(x(globalindex)));
        node1.x0 = p1.x0;
        node1.binary_variables = p1.binary_variables;
        node1.semicont_variables = p1.semicont_variables;
        node1.semibounds = p1.semibounds;
        node1.pid = pid;pid = pid + 1;
        node1.sosgroups = p1.sosgroups;
        node1.sosvariables = p1.sosvariables;
        node1.atmost = p1.atmost;        
                
        node0.lb = p0.lb;
        node0.ub = p0.ub;
        node0.depth = p0.depth;
        node0.lower = p0.lower;
        node0.fixedvariable = index;
        node0.fixdir = 'down';
        node0.TotalIntegerInfeas = TotalIntegerInfeas;
        node0.TotalBinaryInfeas = TotalBinaryInfeas;
        node0.IntInfeas = x(globalindex)-floor(x(globalindex));
        node0.x0 = p0.x0;
        node0.binary_variables = p0.binary_variables;
        node0.semicont_variables = p0.semicont_variables;
        node0.semibounds = p0.semibounds;
        node0.pid = pid;pid = pid + 1;
        node0.sosgroups = p0.sosgroups;
        node0.sosvariables = p0.sosvariables;
        node0.atmost = p0.atmost;        
        
        if ismember(globalindex,p.atmost.variables)
            for j = 1:length(p.atmost.groups)
                xy = p.atmost.groups{j};
                if p.atmost.bounds(j)==1 && any(xy == globalindex)
                    if ~(node0.lb(globalindex)==0 && node0.ub(globalindex)==0)
                        % The variable has been fixed to a non-zero value
                        % Hence, its sister has to be set to 0
                        sisters = xy(xy~=globalindex);
                        for k = sisters
                            if node0.lb(k) > 0 || node0.ub(k) < 0
                                p0_feasible = 0;
                                break
                            else
                                node0.lb(k) = 0;
                                node0.ub(k) = 0;
                            end                            
                        end
                    end
                    if ~(node1.lb(globalindex)==0 && node1.ub(globalindex)==0)
                        sisters = xy(xy~=globalindex);
                        for k = sisters
                            if node1.lb(k) > 0 || node1.ub(k) < 0
                                p1_feasible = 0;
                                break
                            else
                                node1.lb(k) = 0;
                                node1.ub(k) = 0;
                            end                            
                        end
                    end
                end
            end
        end
        
        % Make sure we don't push trivially poor stuff to stack, so reuse
        % pruning code by creating temporary stacks first
        tempstack = stackCreate;
        if p0_feasible            
            tempstack = push(tempstack,node0);
        end
        if p1_feasible
            tempstack = push(tempstack,node1);
        end
        tempstack = prune(tempstack,upper,p.options,solved_nodes,p,allSolutions);
        stack = mergeStack(stack,tempstack);    
    end
        
     if stackLength(stack)>0
        lower = stackLower(stack);
        if can_use_ceil_lower
            lower = ceil(lower);
        end
    end
    
    % Dude, all problems we solve now are infeasible. Start presolving LPs
    % before going all in on the full SDP
    if length(feasibilityHistory) > 5 && all(feasibilityHistory(end-4:end)==0)
        aggresiveprune = 1;
    else
        aggresiveprune = 0;
    end    
    
    % **********************************
    % Get a new node to solve
    % **********************************
    [node,stack] = pull(stack,p.options.bnb.method,x_min,upper);
    if ~isempty(node)
        p = copyNode(p,node);
    end
    
    if isempty(node)
        % There are no nodes left
        if ~isinf(upper)
            gap = 0;
        end
    else
        % We pulled a new node from stack, so there are nodes left
        gap = abs((upper-lower)/(1e-3+abs(upper)+abs(lower)));
    end
    
    if isnan(gap)
        gap = inf;
    end
    
    lowerhist = [lowerhist lower];
    upperhist = [upperhist upper];
    stacksizehist = [stacksizehist stackLength(stack)];
      
    if p.options.bnb.verbose;
        if mod(solved_nodes-1,p.options.print_interval)==0 || isempty(node) || (gap == 0) || (lastUpper-1e-6 > upper)
            if p.options.bnb.plot
                hold off
                subplot(1,3,1);               
                l = plot([lowerhist' upperhist']);set(l,'linewidth',2);
                title('Upper/lower bounds')
                subplot(1,3,2);                
                l = plot(stacksizehist);set(l,'linewidth',2);
                title('Open nodes')
                drawnow
                subplot(1,3,3);   
                hist(getStackLowers(stack),25);
                title('Histogram lower bounds')
                drawnow
            end
            if lastUpper > upper
                fprintf(' %4.0f : %12.3E  %7.2f   %12.3E  %2.0f  %8.1f    %s \n',solved_nodes,upper,100*gap,lower,stackLength(stack),etime(clock,bnbsolvertime),'-> Found improved solution!');
            else
                fprintf(' %4.0f : %12.3E  %7.2f   %12.3E  %2.0f  %8.1f    %s \n',solved_nodes,upper,100*gap,lower,stackLength(stack),etime(clock,bnbsolvertime),yalmiperror(output.problem));
            end
        end
    end
    lastUpper = upper;    
end
if p.options.bnb.verbose;showprogress([num2str2(solved_nodes,3)  ' Finishing.  Cost: ' num2str(upper) ],p.options.bnb.verbose);end
    
% **********************************
%% BRANCH VARIABLE
% **********************************
function [index,whatsplit,globalindex] = branchvariable(x,integer_variables,binary_variables,options,x_min,Weight,p)
all_variables = [integer_variables(:);binary_variables(:)];

if ~isempty(p.sosvariables) && isempty(setdiff(all_variables,p.sosvariables)) & strcmp(options.bnb.branchrule,'sos')
    % All variables are in SOS1 constraints
    for i = 1:length(p.sosgroups)
        dist(i) = (sum(x(p.sosgroups{i}))-max(x(p.sosgroups{i})))/length(p.sosgroups{i});
    end
    % Which SOS to branch on
    [val,index] = max(dist);
    whatsplit = 'sos1';
    globalindex = index;  
end

switch options.bnb.branchrule
    case 'weight'
        interror = abs(x(all_variables)-round(x(all_variables)));
        [val,index] = max(abs(p.weight(all_variables)).*interror);
    case 'first'
        index = min(find(abs(x(all_variables)-round(x(all_variables)))>options.bnb.inttol));
    case 'last'
        index = max(find(abs(x(all_variables)-round(x(all_variables)))>options.bnb.inttol));
    case 'min'
        nint = find(abs(x(all_variables)-round(x(all_variables)))>options.bnb.inttol);
        [val,index] = min(abs(x(nint)));
        index = nint(index);
    case 'max'
        %[val,index] = max((abs(x(all_variables)-round(x(all_variables)))));
        [val,index] = max((1 + min(10,abs(p.c(all_variables)))).*(abs(x(all_variables)-round(x(all_variables)))));
    otherwise
        error('Branch-rule not supported')
end
if index<=length(integer_variables)
    whatsplit = 'integer';
    globalindex = integer_variables(index);
else
    index = index-length(integer_variables);
    whatsplit = 'binary';
    globalindex = binary_variables(index);
end
if isempty(index) | ~isempty(p.semicont_variables)
    for i = 1:length(p.semicont_variables)
        j = p.semicont_variables(i);
        if x(j)>= p.semibounds.lb(i) & x(j)<= p.semibounds.ub(i)
            s(i) = 0;
        elseif x(j)==0
            s(i) = 0;
        else
            s(i) = min([abs(x(j)-0); abs(x(j)-p.semibounds.lb(i));abs(x(j)-p.semibounds.ub(i))]);
        end
    end
    [val2,index2] = max(s);
    if isempty(val)
        whatsplit = 'semi';
        index = index2;
    elseif val2>val
        % index = p.semicont_variables(index);
        whatsplit = 'semi';
        index = index2;
    end
end



% **********************************
% SPLIT PROBLEM
% **********************************
function [p0,p1,variable] = binarysplit(p,x,index,lower,options,sosgroups,sosvariables)
p0 = p;
p1 = p;

variable = p.binary_variables(index);
tf = ~(ismembcYALMIP(p0.binary_variables,variable));
new_binary = p0.binary_variables(tf);

friends = [];
if ~isempty(sosvariables)
    if ismember(variable,sosvariables)
        i = 1;
        while i<=length(sosgroups)            
            if ismember(variable,sosgroups{i})
                friends = setdiff(sosgroups{i},variable);
                break
            else
                i = i + 1;
            end
        end
    end
end

p0.ub(variable)=0;
p0.lb(variable)=0;
if length(friends) == 1
    p0.ub(friends) = 1;
    p0.lb(friends) = 1;
end

p0.lower = lower;
p0.depth = p.depth+1;
p0.binary_variables = new_binary;

p1.ub(variable)=1;
p1.lb(variable)=1;
if length(friends) > 1
    p1.ub(friends)=0;
    p1.lb(friends)=0;
end

p1.binary_variables = new_binary;%p0.binary_variables;%setdiff1D(p1.binary_variables,variable);
%p1.binary_variables = setdiff(p1.binary_variables,friends);
p1.lower = lower;
p1.depth = p.depth+1;

% % *****************************
% % PROCESS MOST PROMISING FIRST
% % (p0 in top of stack)
% % *****************************
if x(variable)>0.5
    pt=p1;
    p1=p0;
    p0=pt;
end

function [p0,p1] = integersplit(p,x,index,lower,options,x_min)

variable = p.integer_variables(index);
current = x(p.integer_variables(index));
lb = floor(current)+1;
ub = floor(current);

% xi<ub
p0 = p;
p0.lower = lower;
p0.depth = p.depth+1;
p0.x0(variable) = ub;
p0.ub(variable)=min(p0.ub(variable),ub);

% xi>lb
p1 = p;
p1.lower = lower;
p1.depth = p.depth+1;
p1.x0(variable) = lb;
p1.lb(variable)=max(p1.lb(variable),lb);

% *****************************
% PROCESS MOST PROMISING FIRST
% *****************************
if lb-current<0.5
    pt=p1;
    p1=p0;
    p0=pt;
end

function [p0,p1] = sos1split(p,x,index,lower,options,x_min)

v = p.sosgroups{index};
n = ceil(length(v)/2);
v1 = v(randperm(length(v),n));
v2 = setdiff(v,v1);

% In first node, set v2 to 0 and v1 to sosgroup
p0 = p;p0.lower = lower;
p0.sosgroups{index} = v1;
p0.ub(v2) = 0;

% In second node, set v1 to 0 and v1 to sosgroup
p1 = p;p1.lower = lower;
p1.sosgroups{index} = v2;
p1.ub(v1) = 0;

function [p0,p1] = semisplit(p,x,index,lower,options,x_min)

variable = p.semicont_variables(index);
current = x(p.semicont_variables(index));

p0 = p;
p0.lower = lower;
p0.depth = p.depth+1;
p0.x0(variable) = 0;
p0.lb(variable)=0;
p0.ub(variable)=0;

p1 = p;
p1.lower = lower;
p1.depth = p.depth+1;
p1.x0(variable) = p.semibounds.lb(index);
p1.lb(variable) = p.semibounds.lb(index);
p1.ub(variable) = p.semibounds.ub(index);

p0.semicont_variables = setdiff(p.semicont_variables,variable);
p1.semicont_variables = setdiff(p.semicont_variables,variable);
p0.semibounds.lb(index)=[];
p0.semibounds.ub(index)=[];
p1.semibounds.lb(index)=[];
p1.semibounds.ub(index)=[];

function s = num2str2(x,d,c);
if nargin==3
    s = num2str(x,c);
else
    s = num2str(x);
end
s = [repmat(' ',1,d-length(s)) s];


function [stack,lower] = prune(stack,upper,options,solved_nodes,p,allSolutions)
% *********************************
% PRUNE STACK W.R.T NEW UPPER BOUND
% *********************************
if stackLength(stack)>0 && ~isinf(upper)
    if length(p.integer_variables) == length(p.c) && all(p.c == fix(p.c)) && nnz(p.Q)==0 && isempty(p.evalMap) && nnz(p.variabletype)==0       
        L = stack.lower;
        tooLarge = find(~isinf(L) & L>=upper-0.999);
    else     
        L = stack.lower;
        tooLarge = find(~isinf(L) & L>=upper*(1-options.bnb.prunetol));
    end
    if ~isempty(tooLarge)        
        stack.nodeCount = stack.nodeCount - length(tooLarge);
        stack.lower(tooLarge) = inf;        
    end
end

% Prune simple linear model w.r.t bound constraints only
if nnz(p.Q) == 0 && isempty(p.evalMap) && nnz(p.variabletype)==0
    tooLarge = [];
    for i = find(~isinf(stack.lower))
        pi = stack.nodes{i};
        neg = find(p.c < 0);
        pos = find(p.c > 0);
        obj = p.f + p.c(pos)'*pi.lb(pos) + p.c(neg)'*pi.ub(neg);
        if obj >= upper
            tooLarge = [tooLarge i];
        end
    end
    if ~isempty(tooLarge)
        stack.nodeCount = stack.nodeCount - length(tooLarge);
        stack.lower(tooLarge) = inf;  
    end    
end

if stack.nodeCount > 0
    lower = min(stack.lower);
else
    lower = upper;
end

function p = adaptivestrategy(p,upper,solved_nodes)
% **********************************'
% SWITCH NODE SELECTION STRATEGY?
% **********************************'
if strcmp(p.options.bnb.method,'depthproject') & (upper<inf)
    p.options.bnb.method = 'project';
end
if strcmp(p.options.bnb.method,'depthbest') & (upper<inf)
    p.options.bnb.method = 'best';
end
if strcmp(p.options.bnb.method,'depthprojection') & (upper<inf)
    p.options.bnb.method = 'projection';
end
if strcmp(p.options.bnb.method,'depthbreadth') & (upper<inf)
    p.options.bnb.method = 'breadth';
end
if strcmp(p.options.bnb.method,'depthest') & (upper<inf)
    p.options.bnb.method = 'est';
end

function res = resids(p,x)
res = [];
if p.K.f>0
    res = -abs(p.F_struc(1:p.K.f,:)*[1;x]);
end
if p.K.l>0
    res = [res;p.F_struc(p.K.f+1:p.K.f+p.K.l,:)*[1;x]];
end
if (length(p.K.q)>1) | p.K.q>0
    top = 1+p.K.f+p.K.l;
    for i = 1:length(p.K.q)
        n = p.K.q(i);
        q = p.F_struc(top:top+n-1,:)*[1;x];top = top+n;
        res = [res;q(1) - norm(q(2:end))];
    end
end
if (length(p.K.s)>1) | p.K.s>0
    top = 1+p.K.f+p.K.l+sum(p.K.q);
    for i = 1:length(p.K.s)
        n = p.K.s(i);
        X = p.F_struc(top:top+n^2-1,:)*[1;x];top = top+n^2;
        X = reshape(X,n,n);
        res = [res;min(eig(X))];
    end
end
res = [res;min([p.ub-x;x-p.lb])];

function [x_min,upper] = initializesolution(p);

x_min = zeros(length(p.c),1);
upper = inf;
if p.options.usex0
    z = p.x0;
    residual = resids(p,z);
    relaxed_feasible = all(residual(1:p.K.f)>=-1e-8) & all(residual(1+p.K.f:end)>=-1e-6);
    if relaxed_feasible & all(z(p.integer_variables)==fix(z(p.integer_variables))) & all(z(p.binary_variables)==fix(z(p.binary_variables)))
        upper = computecost(p.f,p.c,p.Q,z,p);
        x_min = z;
    end
else
    p.x0 = zeros(length(p.c),1);
    x = p.x0;
    z = evaluate_nonlinear(p,x);
    residual = resids(p,z);
    relaxed_feasible = all(residual(1:p.K.f)>=-p.options.bmibnb.eqtol) & all(residual(1+p.K.f:end)>=p.options.bmibnb.pdtol);
    if relaxed_feasible
        upper = computecost(p.f,p.c,p.Q,z,p);
        x_min = x;
    end
end

function p = copyNode(p,node);
p.lb = node.lb;
p.ub = node.ub;
p.depth = node.depth;
p.lower = node.lower;
p.fixedvariable = node.fixedvariable;
p.fixdir = node.fixdir;
p.TotalIntegerInfeas = node.TotalIntegerInfeas;
p.TotalBinaryInfeas = node.TotalBinaryInfeas;
p.IntInfeas =  node.IntInfeas;
p.x0 = node.x0;
p.binary_variables = node.binary_variables;
p.semicont_variables = node.semicont_variables;
p.semibounds = node.semibounds;
p.pid = node.pid;
p.sosgroups = node.sosgroups;
p.sosvariables = node.sosvariables;
p.atmost = node.atmost;

function stack = stackCreate
stack.nodes = {};
stack.lower = [];
stack.nodeCount = 0;

function stack = push(stack,p)
stack.nodes{end + 1} = p;
stack.lower(end + 1) = p.lower;
stack.nodeCount = stack.nodeCount + 1;   

function stack1 = mergeStack(stack1,stack2)
for i = 1:1:length(stack2.nodes)
    if ~isinf(stack2.lower(i))
        stack1.nodes{end + 1} = stack2.nodes{i};
        stack1.lower(end + 1) = stack2.lower(i);
        stack1.nodeCount = stack1.nodeCount + 1;   
    end
end

function stack = compressStack(stack)

used = find(~isinf(stack.lower));
stack.lower = stack.lower(used);
stack.nodes = {stack.nodes{used}};

function [p,stack] = pull(stack,method,x_min,upper);

if stackLength(stack) > 0
    if numel(stack.lower) > 100 && nnz(isinf(stack.lower)) > 0.25*numel(stack.lower)
        stack = compressStack(stack);
    end
    switch method
        case {'depth','depthfirst','depthbreadth','depthproject','depthbest'}
            depths = getStackDepths(stack);
            [i,j] = max(depths);
            % Silly for backward testing compatibility. Stack order has
            % changed, to be able to compare some examples, make sure we
            % traverse the tree in exactly the same was as before on ties
            j = max(find(depths == i));
            p = getStackNode(stack,j);
            stack = removeStackNode(stack,j);
                     
        case 'project'
            error
            [i,j]=min([stack.projection]);
            p=stack(j);
            stack = stack([1:1:j-1 j+1:1:end]);
            
        case 'breadth'
            error
            [i,j]=min([stack.depth]);
            p=stack(j);
            stack = stack([1:1:j-1 j+1:1:end]);
            
        case 'best'                       
            lowers = getStackLowers(stack);
            [i,j] = min(lowers);
            % Silly for backward testing compatibility. Stack order has
            % changed, to be able to compare some examples, make sure we
            % traverse the tree in exactly the same was as before on ties
            j = max(find(lowers == i));
            p = getStackNode(stack,j);
            stack = removeStackNode(stack,j);
            
        otherwise
    end
else
    p = [];
end

function n = stackLength(stack)
n = stack.nodeCount;

function [L,pos,N] = stackLower(stack)
if stack.nodeCount > 0
    if nargout == 1
        L = min(stack.lower);
    elseif nargout == 2
        [L,pos] = min(stack.lower);
    elseif nargout == 3
        [L,pos] = min(stack.lower);
        N = stack.nodes(pos);
    end
else
    L = nan;
end

function D = getStackDepths(stack)
D = -inf(1,length(stack.nodes));
for i = find(~isinf(stack.lower))
    D(i) = stack.nodes{i}.depth;
end

function L = getStackLowers(stack)
L = stack.lower;

function N = getStackNode(stack,j)
N = stack.nodes{j};

function stack = removeStackNode(stack,j)
stack.nodeCount = stack.nodeCount - length(j);
stack.lower(j) = inf;
          
function p = detectSOS(p)
sosgroups = {};
sosvariables = [];
if p.K.f > 0 & ~isempty(p.binary_variables)
    nbin = length(p.binary_variables);
    Aeq = -p.F_struc(1:p.K.f,2:end);
    beq = p.F_struc(1:p.K.f,1);
    notbinary_var_index = setdiff(1:length(p.lb),p.binary_variables);
    only_binary = ~any(Aeq(:,notbinary_var_index),2);
    Aeq_bin = Aeq(find(only_binary),p.binary_variables);
    beq_bin = beq(find(only_binary),:);
    % Detect groups with constraints sum(d_i) == 1
    sosgroups = {};
    for i = 1:size(Aeq_bin,1)
        if beq_bin(i) == 1
            [ix,jx,sx] = find(Aeq_bin(i,:));
            if all(sx == 1)
                sosgroups{end+1} = p.binary_variables(jx);
                sosvariables = [sosvariables p.binary_variables(jx)];
            end
        end
    end
end
p.sosgroups = sosgroups;
p.sosvariables = sosvariables;


function p = simplePresolve(p)
pss=[];
p = propagate_bounds_from_equalities(p);

if p.K.f > 0
    pp = p;
    r = find(p.lb == p.ub);
    pp.F_struc(:,1) = pp.F_struc(:,1) + pp.F_struc(:,r+1)*p.lb(r);
    pp.F_struc(:,r+1)=[];
    pp.lb(r)=[];
    pp.ub(r)=[];
    pp.variabletype(r)=[];
    % FIXME: This is lazy, should update new list
    pp.binary_variables = [];
    pp.integer_variables = [];
    pp = propagate_bounds_from_equalities(pp);
    other = setdiff(1:length(p.lb),r);
    p.lb(other) = pp.lb;
    p.ub(other) = pp.ub;
    p = update_integer_bounds(p);
    redundant = find(~any(pp.F_struc(1:p.K.f,2:end),2));
    if any(p.F_struc(redundant,1)<0)
        p.feasible = 0;
    else
        p.F_struc(redundant,:)=[];
        p.K.f = p.K.f - length(redundant);
    end
end    

function p = extractBounds(p)
if ~isempty(p.F_struc)
    [lb,ub,used_rows_eq,used_rows_lp] = findulb(p.F_struc,p.K);
    if ~isempty(used_rows_lp)
        used_rows_lp = used_rows_lp(~any((p.F_struc(p.K.f + used_rows_lp,1+p.nonlinear)),2));
        if ~isempty(used_rows_lp)
            lower_defined = find(~isinf(lb));
            if ~isempty(lower_defined)
                p.lb(lower_defined) = max(p.lb(lower_defined),lb(lower_defined));
            end
            upper_defined = find(~isinf(ub));
            if ~isempty(upper_defined)
                p.ub(upper_defined) = min(p.ub(upper_defined),ub(upper_defined));
            end
            p.F_struc(p.K.f + used_rows_lp,:)=[];
            p.K.l = p.K.l - length(used_rows_lp);
        end
    end
    
    if ~isempty(used_rows_eq)
        used_rows_eq = used_rows_eq(~any(full(p.F_struc(used_rows_eq,1+p.nonlinear)),2));
        if ~isempty(used_rows_eq)
            lower_defined = find(~isinf(lb));
            if ~isempty(lower_defined)
                p.lb(lower_defined) = max(p.lb(lower_defined),lb(lower_defined));
            end
            upper_defined = find(~isinf(ub));
            if ~isempty(upper_defined)
                p.ub(upper_defined) = min(p.ub(upper_defined),ub(upper_defined));
            end
            p.F_struc(used_rows_eq,:)=[];
            p.K.f = p.K.f - length(used_rows_eq);
        end
    end
end

function p = detect3x3SymmetryGroups(p)

good = zeros(1,length(p.c));
good(p.integer_variables) = 1;
good( (p.lb ~= -1) | (p.ub ~=0)) = 0;
if any(good)
    groups = {};
    for j = 1:length(p.K.s)
        if p.K.s(j) >= 3
        n = 3;
        X = spalloc(p.K.s(j),p.K.s(j),p.K.s(j)^2);
        X(1:n,1:n) = 1;
        index0 = find(X);
        index = index0;
        corner = 0;        
        for block = 1:p.K.s(j)-n
            dataBlock = p.semidefinite{j}.F_struc(index,:);
            used = find(any(dataBlock,1));
            dataBlock = dataBlock(:,used);             
            if used(1) == 0;dataBlock = [zeros(size(dataBlock,1),1) dataBlock];end
            v = used;v = v(v>1)-1;
            if all(good(v))
                if isempty(groups)
                    groups{1}.dataBlock = dataBlock;
                    groups{1}.variables{1} = used;
                else
                    found = 0;
                    for k = 1:length(groups)
                        if isequal(length(used),length(groups{1}.variables{1}))
                            if isequal(groups{1}.dataBlock,dataBlock)
                                found = 1;
                                groups{1}.variables{end+1} = used;
                            else
                                % TODO: Look for simple scaled versions
                                % s = groups{1}.dataBlock(:)\dataBlock(:);
                                % norm(s*groups{1}.dataBlock-dataBlock,inf)
                            end
                        end
                    end
                    if ~found
                        groups{end+1}.dataBlock = dataBlock;
                        groups{end}.variables{1} = used;
                    end
                end
            end
            index = index + p.K.s(j)+1;
        end
        end
    end
    for i = 1:length(groups)
        if length(groups{i}.variables) > 1
            keep(i) = 1;
        else
            keep(i) = 0;
        end
    end
    if length(groups) > 0
        groups = {groups{find(keep)}};
        if length(groups) > 0
            for i = 1:length(groups)
                for j = 1:length(groups{i}.variables);
                    v = groups{i}.variables{j};
                    v = v(v>1)-1;
                    groups{i}.variables{j} = v;
                end
            end
        end
    end
else
    groups = {};
end
p.sdpsymmetry = groups;
    
    
function p_lp = add3x3sdpsymmetrycut(p,p_lp,x) 

for j = 1:length(p.sdpsymmetry)    
    excludes = [];
    n = sqrt(size(p.sdpsymmetry{j}.dataBlock,1));
    for i = 1:length(p.sdpsymmetry{j}.variables)        
        if min(eig(reshape(p.sdpsymmetry{j}.dataBlock*[1;x(p.sdpsymmetry{j}.variables{i})],n,n))) < -abs(p_lp.options.bnb.feastol)
            excludes = [excludes x(p.sdpsymmetry{j}.variables{i})];
        end
    end
    if ~isempty(excludes)
        newF = [];
        infeasible_combinations = unique(excludes','rows')';
        for k = 1:size(infeasible_combinations,2)
            % Local cut for reduced set
            [b,atemp] = exclusionCut(infeasible_combinations(:,k),-1);
            % Add that cut for every variable groups
            for s = 1:length(p.sdpsymmetry{j}.variables)
                a = spalloc(1,length(p_lp.c),1);
                a(p.sdpsymmetry{j}.variables{s}) = atemp;
                if all(sum(abs(p_lp.F_struc - [b a]),2)<=1e-12)
                    newF = [newF;b a];
                end
            end
        end
        p_lp = addLinearCut(p_lp,newF);
    end
end
    
    
function [p,p_lp] = addSymmetryCuts(p,p_lp)

for j = 1:length(p.sdpsymmetry)    
    if length(p.sdpsymmetry{j}.variables{1}) <= 4
        % We can enumerate easily infeasibles
        excludes = [];
        n = sqrt(size(p.sdpsymmetry{j}.dataBlock,1));
        combs = -dec2decbin(0:2^length(p.sdpsymmetry{j}.variables{1})-1,length(p.sdpsymmetry{j}.variables{1}))';
        for i = 1:size(combs,2)
            if min(eig(reshape(p.sdpsymmetry{j}.dataBlock*[1;combs(:,i)],n,n))) < -abs(p_lp.options.bnb.feastol)
                excludes = [excludes combs(:,i)];
            end
        end
        if ~isempty(excludes)
            newF = [];
            infeasible_combinations = unique(excludes','rows')';
            for k = 1:size(infeasible_combinations,2)
                % Local cut for reduced set
                [b,atemp] = exclusionCut(infeasible_combinations(:,k),-1);
                % Add that cut for every variable groups
                for s = 1:length(p.sdpsymmetry{j}.variables)
                    a = spalloc(1,length(p_lp.c),1);
                    a(p.sdpsymmetry{j}.variables{s}) = atemp;                    
                    newF = [newF;b a];                  
                end
            end
            p_lp = addLinearCut(p_lp,newF);
            % Delete, won't need these in the future
            p.sdpsymmetry{j}.dataBlock = [];
            p.sdpsymmetry{j}.variables = [];
        end
    end
end

function p = addLinearCut(p,row);
if sum(p.K.s) == 0 && sum(p.K.q)==0
    % Just append in the end
    p.F_struc = [p.F_struc;row];
else
    % Insert before conics
    p.F_struc = [p.F_struc(1:(p.K.f+p.K.l),:);row;p.F_struc(1+p.K.f+p.K.l:end,:)];
end
p.K.l = p.K.l + size(row,1);