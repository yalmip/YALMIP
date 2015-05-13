function output = callqpas(interfacedata)

% Retrieve needed data
options = interfacedata.options;
model = yalmip2quadprog(interfacedata);

if size(model.Aeq,1)==0
    model.Aeq = [];
    model.beq = [];
end

if options.savedebug
    save debugfile model
end

if options.showprogress;showprogress(['Calling ' interfacedata.solver.tag],options.showprogress);end
solvertime = tic;
[x,flag,lm] = qpas(model.Q, model.c, model.A, model.b, model.Aeq, model.beq, model.lb, model.ub, options.verbose);
solvertime = toc(solvertime);

% Internal format for duals
if ~isempty(lm)
    D_struc = [lm.equality;lm.inequality];
else
    D_struc = [];
end

if isempty(flag)
    problem = 9;
else
    % Check, currently not exhaustive...
    switch flag
        case 0
            problem = 0;
        case 2
            problem = 1;
        otherwise
            problem = 1;
    end
end

infostr = yalmiperror(problem,interfacedata.solver.tag);

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
output = createOutputStructure(x(:),D_struc,[],problem,infostr,solverinput,solveroutput,solvertime);