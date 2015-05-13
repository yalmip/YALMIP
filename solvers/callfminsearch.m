function output = callfminsearch(model)

% Retrieve needed data
options = model.options;
c       = model.c;
x0      = model.x0;
Q       = model.Q;
lb      = model.lb;
ub      = model.ub;
K       = model.K;
monomtable = model.monomtable;

nobounds = isempty(lb) | (~isempty(lb) & all(isinf(lb)));
nobounds = nobounds & (isempty(ub) | (~isempty(ub) & all(isinf(ub))));
if any(K.l) | any(K.q) | any(K.s) | ~nobounds
    output = createoutput(zeros(length(c),1),[],[],-4,'FMINSEARCH',[],[],0);
    return
end

if isempty(model.evaluation_scheme)
    model = build_recursive_scheme(model);
end

switch options.verbose
    case 0
        options.fmincon.Display = 'off';
    case 1
        options.fmincon.Display = 'final';
    otherwise
        options.fmincon.Display = 'iter';
end

% Do some pre-calc to be used in calls from fmincon
nonlinearindicies = union(find(model.variabletype~=0),model.evalVariables);
linearindicies    = setdiff(find(model.variabletype==0),nonlinearindicies);
model.nonlinearindicies = nonlinearindicies;
model.linearindicies    = linearindicies;

model.Anonlinineq = [];
model.bnonlinineq = [];
model.Anonlineq = [];
model.bnonlineq = [];

% Extract linear and nonlinear equality constraints
Aeq = [];
beq = [];
A = [];
b = [];

% This helps with robustness in bnb in some cases
x0candidate = zeros(length(c),1);

if isempty(x0)
    x0 = x0candidate(linearindicies);
else
    x0 = x0(linearindicies);
end

if options.savedebug
    ops = options.fminsearch;
    save fminsearch model x0 ops
end

showprogress('Calling FMINSEARCH',options.showprogress);

% Precalc for the callbacks
model = setup_fmincon_params(model);
if (model.SimpleQuadraticObjective | model.SimpleNonlinearObjective) & isempty(model.evalMap)
    options.fmincon.GradObj = 'on';    
end

fun = @(x,m)fmincon_fun(x,m);

solvertime = tic;
[xout,fmin,flag,output] = fminsearch(@(x)fun(x,model),x0,options.fminsearch);
solvertime = toc(solvertime);

if isempty(nonlinearindicies)
    x = xout(:);
else
    x = zeros(length(c),1);
    for i = 1:length(linearindicies)
        x(linearindicies(i)) = xout(i);
    end
    x = x(1:length(c));
end

problem = 0;

% Internal format for duals
D_struc = [];

% Check, currently not exhaustive...
if flag==0
    problem = 3;
else
    problem = 0;
end

% Save all data sent to solver?
if options.savesolverinput   
    solverinput.fun = 'fmincon_fun';
    solverinput.x0 = x0;
    solverinput.ops = options.fminsearch;
    solverinput.param = model;
else
    solverinput = [];
end

% Save all data from the solver?
if options.savesolveroutput
    solveroutput.x = x;
    solveroutput.fmin = fmin;
    solveroutput.flag = flag;
    solveroutput.output=output;
    solveroutput.lambda=lambda;
else
    solveroutput = [];
end

% Standard interface
output = createoutput(x,D_struc,[],problem,'FMINSEARCH',solverinput,solveroutput,solvertime);