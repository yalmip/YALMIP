function output = callknitro_qp(interfacedata)

model = yalmip2knitromiqp(interfacedata);

showprogress('Calling KNITRO',interfacedata.options.showprogress);
solvertime = tic;
[x,fval,exitflag,output,lambda] = knitro_qp(model.H,model.f,model.A,model.b,model.Aeq,model.beq,model.lb,model.ub,model.x0,model.ops);
solvertime = toc(solvertime);

% Internal format for duals
D_struc = [lambda.eqlin;lambda.ineqlin];

% Check, currently not exhaustive...
problem = 0;
switch exitflag
    case {0,-101}
        problem = 0;
    case {-200,-204,-205,-515}
        problem = 1;
    case {-101,-300}
        problem = 2;
    case {-202,-400,-401,-410}
        problem = 3;
    otherwise
        problem = 11;
end

% Save all data sent to solver?
if interfacedata.options.savesolverinput
    solverinput.model = model;   
else
    solverinput = [];
end

% Save all data from the solver?
if interfacedata.options.savesolveroutput
    solveroutput.x = x;
    solveroutput.fval = fval;
    solveroutput.exitflag = exitflag;
    solveroutput.output=output;
    solveroutput.lambda=lambda;
else
    solveroutput = [];
end

% Standard interface
output = createOutputStructure(x,D_struc,[],problem,interfacedata.solver.tag,solverinput,solveroutput,solvertime);