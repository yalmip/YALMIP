function output = call_cplexibm_qcmiqp(interfacedata)

% Author Johan Löfberg

% This is a gateway to all CPLEX interfaces
% Call LP/QP solver if sufficient
% Bug in 12.6 regarding display options. Hence, we turn it on and then run
% silently by redirecting output instead
interfacedata.options.cplex.Display = 'on';
if isempty(interfacedata.K.q) | interfacedata.K.q(1)==0
    output = call_cplexibm_miqp(interfacedata);
    return
end

options = interfacedata.options;
n_original = length(interfacedata.c);
[model,nonlinearremain] = yalmip2cplex(interfacedata);
if nonlinearremain
    error('Nonlinear monomials remain when calling CPLEX. If you are using OPTIMIZER, ensure your model really is solvable by CPLEX for fixed parameters. If you still think so, please report this and ask for a feature improvement.');
end
if options.savedebug
    save cplexdebug model
end

% Call mex-interface
showprogress('Calling CPLEX',options.showprogress);
if isempty(model.integer_variables) & isempty(model.binary_variables) & isempty(model.semicont_variables) & isempty(model.K.sos.type)
    if options.verbose
        solvertime = tic;
        [x,fval,exitflag,output] = cplexqcp(model.H, model.f, model.Aineq,model.bineq,model.Aeq,model.beq,model.Li,model.Qi,model.ri,model.lb,model.ub,model.x0,model.options);
        solvertime = toc(solvertime);
    else
        solvertime = tic;
        evalc('[x,fval,exitflag,output] = cplexqcp(model.H, model.f, model.Aineq,model.bineq,model.Aeq,model.beq,model.Li,model.Qi,model.ri,model.lb,model.ub,model.x0,model.options);');
        solvertime = toc(solvertime);
    end
else
    if options.verbose
        solvertime = tic;
        [x,fval,exitflag,output] = cplexmiqcp(model.H, model.f, model.Aineq,model.bineq,model.Aeq,model.beq,model.Li,model.Qi,model.ri,model.K.sos.type,model.K.sos.variables,model.K.sos.weight,model.lb,model.ub,model.ctype,model.x0,model.options);
        solvertime = toc(solvertime);
    else
        solvertime = tic;
        evalc('[x,fval,exitflag,output] = cplexmiqcp(model.H, model.f, model.Aineq,model.bineq,model.Aeq,model.beq,model.Li,model.Qi,model.ri,model.K.sos.type,model.K.sos.variables,model.K.sos.weight,model.lb,model.ub,model.ctype,model.x0,model.options);');
        solvertime = toc(solvertime);
    end
end

if length(x) == length(model.f)
    if ~isempty(model.NegativeSemiVar)
        x(model.NegativeSemiVar) = -x(model.NegativeSemiVar);
    end
end

if isempty(x)
    x = zeros(n_original,1);
else
    x = x(1:n_original);
end

problem = 0;
D_struc = [];

% Check, currently not exhaustive...
switch output.cplexstatus
    case {1,101,102}
        problem = 0;
    case {3,103,106}
        problem = 1; % Infeasible
    case {2,20,21,118}
        problem = 2; % Unbounded
    case 4
        problem = 1;
    case {10,11,104,105,107,108,111,112}
        problem = 3; % Iteration/time
    case {5,6,109,110}
        problem = 4; % Numerics
    case 119
        problem = 15;
    otherwise
        problem = -1;
end

infostr = yalmiperror(problem,'CPLEX-IBM');

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
    solveroutput.output=output;
else
    solveroutput = [];
end

% Standard interface
output = createOutputStructure(x,D_struc,[],problem,infostr,solverinput,solveroutput,solvertime);