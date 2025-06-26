function output = callpiqp(interfacedata)

options = interfacedata.options;
model = yalmip2piqp(interfacedata);

if options.savedebug
    save debugfile model
end

if options.showprogress;showprogress(['Calling ' interfacedata.solver.tag],options.showprogress);end

% Define verbose option
options.piqp.verbose = options.verbose;

% Determine whether problem is sparse
is_sparse = issparse(model.P) || issparse(model.A) || issparse(model.G);

% Solve with PIQP
if is_sparse
    PIQPSolver = piqp('sparse');
else
    PIQPSolver = piqp('dense');
end
PIQPSolver.setup(model.P, model.c, ...
                 model.A, model.b, ...
                 model.G, model.h, ...
                 model.x_lb, model.x_ub, ...
                 options.piqp);
results = PIQPSolver.solve();

switch results.info.status_val
    case 1
        problem = 0;
    case -1
        problem = 3;
    case -2
        problem = 1;
    case -3
        problem = 1;
    case -8
        problem = 4;
    case -9
        problem = 7;
    case -10
        problem = 7;
    otherwise
        problem = -1;
end

% Solver time
solvertime = results.info.run_time;

% Standard interface
Primal      = results.x(:);
Dual        = [results.y; results.z];
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
