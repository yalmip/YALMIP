function output = call_xpressfico_milp(interfacedata)

% Retrieve needed data
options = interfacedata.options;
model = yalmip2xpress(interfacedata);

if options.savedebug
    save xpressdebug model
end

solvertime = tic;
if isempty(model.extra.integer_variables) & isempty(model.extra.binary_variables) & isempty(model.extra.semicont_variables) & isempty(model.sos)
    if options.verbose
        [x,fval,exitflag,output,lambda] = xprslp(model.f,model.A,model.b,model.rtype,model.lb,model.ub,model.ops.xpress);
    else
        evalc('[x,fval,exitflag,output,lambda] = xprslp(model.f,model.A,model.b,model.rtype,model.lb,model.ub,model.ops.xpress);');
    end
else
    if options.verbose
        [x,fval,exitflag,output] = xprsmip(model.f,model.A,model.b,model.rtype,model.ctype,model.clim,model.sos,model.lb,model.ub,[],model.ops.xpress);
    else
        evalc('[x,fval,exitflag,output] = xprsmip(model.f,model.A,model.b,model.rtype,model.ctype,model.clim,model.sos,model.lb,model.ub,[],model.ops.xpress);        ');
    end
    lambda = [];
end
solvertime = toc(solvertime);

if ~isempty(lambda)
    D_struc = [lambda.lin];
else
    D_struc = [];
end

if length(x) == length(model.f)
    if ~isempty(model.extra.NegatedSemiVar)
        x(model.extra.NegatedSemiVar) = -x(model.extra.NegatedSemiVar);
    end
end

% Check, currently not exhaustive...
switch exitflag
    case {1}
        problem = 0;
    case {-2}
        problem = 1; % Infeasible
    case {0,-4,-5}
        problem = 3;
    case -8
        problem = -11;
    otherwise
        problem = -1;
end
infostr = yalmiperror(problem,'FICO Xpress');

% Save all data sent to solver?
if options.savesolverinput
    solverinput = model;
else
    solverinput = [];
end

% Save all data from the solver?
if options.savesolveroutput
    solveroutput.x = x;
    solveroutput.fval = fval;
    solveroutput.exitflag = exitflag;
    solveroutput.output = output;
    solveroutput.lambda = lambda;
else
    solveroutput = [];
end

if isempty(x)
    x = zeros(length(model.f),1);
end

% Standard interface
output = createOutputStructure(x(:),D_struc,[],problem,infostr,solverinput,solveroutput,solvertime);