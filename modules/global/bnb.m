function output = bnb(p)
%BNB General branch-and-bound scheme for (primarily) conic programs
%
% BNB applies a branch-and-bound scheme to solve mixed-integer
% convex programs, in particular linear mixed-integer semidefinite
% programs
%
% BNB is never called by the user directly, but is called by
% YALMIP from OPTIMIZE, by choosing the solver tag 'bnb' in sdpsettings.
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

bnbsolvertime = clock;

% Remove options if none has been changed
% Improves performance when calling solver many times
p.options = pruneOptions(p.options);

% We might have a GP : pre-calc
p.nonlinear = find(~(sum(p.monomtable~=0,2)==1 & sum(p.monomtable,2)==1));
p.nonlinear = union(p.nonlinear,p.evalVariables);

% This field is only used in bmibnb, which uses the same sub-functions as
% bnb
p.high_monom_model = [];

% Define infinite bounds, assume feasible
if isempty(p.ub)
    p.ub = repmat(inf,length(p.c),1);
end
if isempty(p.lb)
    p.lb = repmat(-inf,length(p.c),1);
end
p.feasible = 1;

% Various classifications and detection of special stuff
p.rowtype = zeros(1,p.K.f+p.K.l);
p.cliques.table = spalloc(length(p.c),0,0);
p.cliques.hash  = spalloc(length(p.c),0,0);
p.knapsack.a = {};
p.knapsack.b = {};
p.knapsack.variables = {};
p.knapsack.type = [];
p.knapsack.hash = [];
p.hash = randn(length(p.c),1);
p.isbinary  = ismember(1:length(p.c),p.binary_variables);
p.isinteger = ismember(1:length(p.c),p.integer_variables);
p.isintegral = p.isbinary | p.isinteger;

% Could be some nonlinear terms 
% (although these problems are recommended to be solved using BMIBNB
p = compile_nonlinear_table(p);

% Extract bounds from model
p = presolve_bounds_from_modelbounds(p,1);
p = bounds_from_cones_to_lp(p);

% Data loaded as integer, but trivially a binary?
p = presolve_binaryinteger(p);

% Wrong encoding in some historical files with integer -1 -> 0
[p,negated_binary] = detect_negated_binary(p);
p.integral_variables = union(p.binary_variables,p.integer_variables);
p.noninteger_variables = setdiff(1:length(p.c),[p.integral_variables p.semicont_variables]);

% We keep track of explicit cardinality bounds for some propagations
p.binarycardinality.up = length(p.binary_variables);
p.binarycardinality.down = 0;

% Say hello
if p.options.verbose
    disp('* Starting YALMIP integer branch & bound.');
    disp(['* Lower solver   : ' p.solver.lower.tag]);
    disp(['* Upper solver   : ' p.options.bnb.uppersolver]);
    disp(['* Max time       : ' num2str(p.options.bnb.maxtime)]);
    disp(['* Max iterations : ' num2str(p.options.bnb.maxiter)]);      
    nB = length(p.binary_variables);
    nI = length(p.integer_variables);
    nC = length(p.c) - nB - nI;
    disp(['* Inital model: ' num2str(nB) ' binaries, ' num2str(nI) ' integers, ' num2str(nC) ' continuous.']);    
    disp(['*               ' num2str(p.K.f+p.K.l) ' linear rows, ' num2str(sum(p.K.q)+sum(p.K.s.^2)) ' conic rows']);    
    disp(['* Starting presolve and model analysis']);
end

% Clean up model and extract basics
% We remove constraints, but never remove any variables. Instead, we simply
% simply zero out the columns and move to constant
p = classifyrows(p);
p = presolve_empty_rows(p);
p = update_integer_bounds(p);
p = presolve_euclideanreduction(p); 
p = presolve_implied_from_SDP(p);
p = presolve_remove_repeatedrows(p);
p = presolve_remove_dominatedinequalities(p);
p = presolve_remove_equalitydominateinequality(p);
p = presolve_impliedintegers_from_equalities(p);
p = presolve_strengthen_coefficients(p);
p = presolve_dualreductions(p);
p = update_semicont_bounds(p);
p = smashFixed(p);
p = presolve_empty_rows(p);
p = detectRedundantInfeasibleSDPRows(p);
p = presolve_empty_rows(p);
p = presolve_implied_binaryproduct(p);
p = presolve_downforce(p);
p = presolve_upforce(p);
% These are reused from BMIBNB and involves nonlinear models
p = propagate_bounds_from_monomials(p);
p = propagate_bounds_from_equalities(p);
p = propagate_bounds_from_monomials(p);
p = propagate_bounds_from_equalities(p);
p = propagate_bounds_from_monomials(p);
p = propagate_bounds_from_equalities(p);

% Some tricks and strategies are performed for certain objectives
p = detect_special_objectives(p);

if ~p.feasible || any(p.lb > p.ub+1e-6)
    output = createOutputStructure([],[],[],1,yalmiperror(1,'BNB'),[],[],0);    
    return
end

%% ************************************************************************
% Display logics
% 0 : Silent
% 1 : Display branching
% 2 : More YALMIP info
% 3 : Display node solver prints
% *************************************************************************
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

% Figure out the sos weights if any (experimental)
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


% START BRANCHING
setuptime = etime(clock,bnbsolvertime);
bnbsolvertime = clock;
[x_min,solved_nodes,lower,upper,profile,diagnostics] = bnb_branch_and_bound(p);
bnbsolvertime =  etime(clock,bnbsolvertime);
output.solvertime = setuptime + bnbsolvertime;

% Map back in weird models with negated binary
if ~isempty(x_min) && ~isempty(negated_binary)
    x_min(negated_binary) = -x_min(negated_binary);
end

% CREATE SOLUTION
if diagnostics == -4
    output.problem = -4;
elseif diagnostics == 9
    output.problem = 9;
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
    output.solveroutput.localsolver = p.solver.lower.tag;
    output.solveroutput.branchingtime = bnbsolvertime;
    output.solveroutput.solved_nodes = solved_nodes;
    output.solveroutput.lower = lower;
    output.solveroutput.upper = upper;
else
    output.solveroutput =[];
end