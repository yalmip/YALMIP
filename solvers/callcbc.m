function output = callcbc(interfacedata)

% Standard input interface
model = yalmip2cbc(interfacedata);

options = interfacedata.options;
model.opts = options.cbc;
model.opts.display = options.verbose;

if options.savedebug
    save cbcdebug model
end

showprogress('Calling CBC',options.showprogress);
solvertime = tic;
if nnz(model.H)==0
    [x,fval,exitflag,nodes,xc] = cbc([],model.f, model.A, model.rl, model.ru, model.lb, model.ub, model.xtype,model.sos,model.x0,model.opts);
else
    [x,fval,exitflag,nodes,xc] = cbc(model.H,model.f, model.A, model.rl, model.ru, model.lb, model.ub, model.xtype,model.sos,model.x0,model.opts);
end
solvertime = toc(solvertime);

% No duals
D_struc = [];
if isempty(x)
    x = zeros(length(c),1);
end

% Check, currently not exhaustive...
problem = 0;
switch exitflag
    case {0,2,6}
        problem = 0;
    case {1,8}
        problem = 1;
    case {3,4}
        problem = 3;
    case 5
        problem = 16;
    case 7
        problem = 2;
    otherwise
        problem = -1;
end

% Save all data sent to solver?
if options.savesolverinput
    solverinput.model = model;
else
    solverinput = [];
end

% Save all data from the solver?
if options.savesolveroutput
    solveroutput.x = x;
    solveroutput.fval = fval;
    solveroutput.exitflag = exitflag;
    solveroutput.nodes = nodes;
    solveroutput.xc = xc;
else
    solveroutput = [];
end

% Standard interface
output = createOutputStructure(x(:),[],[],problem,interfacedata.solver.tag,solverinput,solveroutput,solvertime);