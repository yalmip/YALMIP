function output = callqpoases(interfacedata)

options = interfacedata.options;
model = yalmip2quadprog(interfacedata);

if options.savedebug
    save debugfile model
end

if options.showprogress;showprogress(['Calling ' interfacedata.solver.tag],options.showprogress);end
solvertime = clock;
solution = callsolver(model,options);
if interfacedata.getsolvertime solvertime = etime(clock,solvertime);else solvertime = 0;end

% Save all data sent to solver?
if ~options.savesolverinput
    model = [];
end

% Save all data from the solver?
if ~options.savesolveroutput
    solveroutput = [];
end

switch solution.exitflag
    case 0
        solution.problem = 0;
    case 1
        solution.problem = 3;
    case -2
        solution.problem = 1;
    case -3
        solution.problem = 2;
    case -1
        solution.problem = 11;
    otherwise
        solution.problem = -1;
end

% Standard interface
output.Primal      = solution.x(:);
output.Dual        = solution.lambda(:);
output.Slack       = [];
output.problem     = solution.problem;
output.infostr     = yalmiperror(solution.problem,interfacedata.solver.tag);
output.solvertime  = solvertime;
if ~options.savesolverinput
    output.solverinput = [];
else
    output.solverinput = model;
end
if ~options.savesolveroutput
    output.solveroutput = [];
else
    output.solveroutput = solution;
end

function solveroutput = callsolver(model,options)
x = [];
fval = [];
exitflag = [];
iter = [];
lambda = [];
lbA = [model.beq;-inf(length(model.b),1)];
ubA = [model.beq;model.b];
A = [model.Aeq;model.A];
options.qpoases.printLevel = -1;
[x,fval,exitflag,iter,lambda] = qpOASES(model.Q, model.c, A, model.lb,model.ub,lbA,ubA,options.qpoases);
solveroutput.x = x;
solveroutput.fval = fval;
solveroutput.exitflag = exitflag;
solveroutput.iter = iter;
if ~isempty(lambda)
solveroutput.lambda = -lambda(1+length(model.lb):end);
else
    solveroutput.lambda=[];
end
