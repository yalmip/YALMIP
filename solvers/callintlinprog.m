function output = callintlinprog(interfacedata)

% Standard input interface
options = interfacedata.options;
model = yalmip2intlinprog(interfacedata);

if options.savedebug    
    save intlinprogdebug model ops
end
 
solvertime = tic;
showprogress('Calling INTLINPROG',options.showprogress);
[x,fval,exitflag,output] = intlinprog(model.c, model.intcon, model.A, model.b, model.Aeq, model.beq, model.lb, model.ub,model.ops);
solvertime = toc(solvertime);
problem = 0;

if isempty(x)
    x = zeros(length(model.c),1);
end

% Check, currently not exhaustive...
switch exitflag
    case 1
        problem = 0;
    case -2
        problem = 1;
    case -3
        problem = 2;       
    otherwise
        problem = 9;        
end
infostr = yalmiperror(problem,'INTLINPROG');       

% Save all data sent to solver?
if options.savesolverinput
    solverinput.model = model;    
    solverinput.options = ops;
else
    solverinput = [];
end

% Save all data from the solver?
if options.savesolveroutput
    solveroutput.x = x;
    solveroutput.fval = fval;
    solveroutput.exitflag = exitflag;
    solveroutput.output=output;
else
    solveroutput = [];
end

% Standard interface 
output = createOutputStructure(x(:),[],[],problem,infostr,solverinput,solveroutput,solvertime);