function output = callqpip(interfacedata)

% Retrieve needed data
options = interfacedata.options;
model = yalmip2quadprog(interfacedata);

if options.savedebug
    model.ops = options.qpip;
    save debugfile model
end

if options.showprogress;showprogress(['Calling ' interfacedata.solver.tag],options.showprogress);end
solvertime = tic;
[x,flag,lm] = qpip(model.Q, model.c, model.A, model.b, model.Aeq, model.beq, model.lb, model.ub, options.verbose,options.qpip.mu,options.qpip.method);
solvertime = toc(solvertime);

% Internal format for duals
if ~isempty(lm)
    D_struc = [lm.equality;lm.inequality];
else
    D_struc = [];
end

% Check, currently not exhaustive...
switch flag
    case 0
        problem = 0;
    case 2
        problem = 1;
    otherwise
        problem = 1;
end

% Save all data sent to solver?
if options.savesolverinput
    solverinput = model;
else
    solverinput = [];
end

% Save all data from the solver?
if options.savesolveroutput
    solveroutput.x = x;
    solveroutput.flag = flag;
    solveroutput.lm = lm;
else
    solveroutput = [];
end

% Standard interface 
output = createOutputStructure(x(:),D_struc,[],problem,interfacedata.solver.tag,solverinput,solveroutput,solvertime);