function output = callosqp(interfacedata)

options = interfacedata.options;
model = yalmip2quadprog(interfacedata);

if options.savedebug
    save debugfile model
end

if options.showprogress;showprogress(['Calling ' interfacedata.solver.tag],options.showprogress);end

% Define QP
n_var = length(model.c);
P = model.Q;
q = model.c;
eye_n = eye(n_var);
A = [model.Aeq;model.A; eye_n];
l = full([model.beq; -inf(length(model.b),1); model.lb]);
u = full([model.beq; model.b; model.ub]);

% Define verbose option
options.osqp.verbose = options.verbose;

% Solve with OSQP
OSQPSolver = osqp;
OSQPSolver.setup(P, q, A, l, u, options.osqp);
results = OSQPSolver.solve();

switch results.info.status_val
    case 1
        problem = 0;
    case -2
        problem = 3;
    case -3
        problem = 1;
    case -4
        problem = 2;
    case -5
        problem = 16;
    otherwise
        problem = -10;
end

% Solver time
solvertime = results.info.run_time;

% Standard interface
Primal      = results.x(:);
Dual        = results.y(1:end-n_var);
infostr     = yalmiperror(problem,interfacedata.solver.tag);
if ~options.savesolverinput
    solverinput = [];
else
    solverinput = model;
end
if ~options.savesolveroutput
    solveroutput = [];
else
    solveroutput = results;
end

% Standard interface
output = createOutputStructure(Primal,Dual,[],problem,infostr,solverinput,solveroutput,solvertime);
