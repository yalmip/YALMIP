function output = call_xpressfico_qcmip(interfacedata)

if nnz(interfacedata.Q) == 0 && nnz(interfacedata.K.q)==0
    output = call_xpressfico_milp(interfacedata);
    return
elseif nnz(interfacedata.K.q)==0
    output = call_xpressfico_miqp(interfacedata);
    return    
end

% Retrieve needed data
options = interfacedata.options;
n_original = length(interfacedata.c);
model = yalmip2xpress(interfacedata);

if options.savedebug
    save xpressdebug model
end

solvertime = tic;
if isempty(model.extra.integer_variables) & isempty(model.extra.binary_variables) & isempty(model.extra.semicont_variables) & isempty(model.sos)
    if options.verbose
        [x,fval,exitflag,output,lambda] = xprsqcqp(model.H,model.f,model.A,model.Q,model.b,model.rtype,model.lb,model.ub,model.ops);
    else
        evalc('[x,fval,exitflag,output,lambda] = xprsqcqp(model.H,model.f,model.A,model.Q,model.b,model.rtype,model.lb,model.ub,model.ops);');
    end
else
    if options.verbose
        [x,fval,exitflag,output] = xprsqcqp(model.H,model.f,model.A,model.Q,model.b,model.rtype,model.ctype,model.clim,model.sos,model.lb,model.ub,[],model.ops);
    else
        evalc('[x,fval,exitflag,output] = xprsqcqp(model.H,model.f,model.A,model.Q,model.b,model.rtype,model.ctype,model.clim,model.sos,model.lb,model.ub,[],model.ops);        ');
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
    x = x(1:n_original);
else
    x = zeros(n_original,1);
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

% Save all data sent to solver?
if options.savesolverinput
    solverinput.f = f;
    solverinput.A = A;
    solverinput.b = b;
    solverinput.ctype = ctype;
    solverinput.rtype = rtype;
    solverinput.lb = lb;
    solverinput.ub = ub;
    solverinput.x0 = [];
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
output = createOutputStructure(x(:),D_struc,[],problem,interfacedata.solver.tag,solverinput,solveroutput,solvertime);