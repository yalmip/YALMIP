function output = callknitro(model)

% Standard NONLINEAR setup
model = yalmip2nonlinearsolver(model);

model.options.knitro.GradObj = 'on';  
model.options.knitro.GradConstr = 'on';
model.options.knitro.JacobPattern = jacobiansparsityfromnonlinear(model,0);

% If quadratic objective and no nonlinear constraints, we can supply an
% Hessian of the Lagrangian
usedinObjective = find(model.c | any(model.Q,2));
if ~any(model.variabletype(usedinObjective)) & any(model.Q)
    if  length(model.bnonlinineq)==0 & length(model.bnonlineq)==0
        H = model.Q(:,model.linearindicies);
        H = H(model.linearindicies,:);
        model.options.knitro.Hessian = 'user-supplied';
        model.options.knitro.HessPattern = sparse(H | H);
        model.options.knitro.HessFcn = @(x,l) 2*H;
    end
end

global latest_xevaled
global latest_x_xevaled
latest_xevaled = [];
latest_x_xevaled = [];

showprogress('Calling KNITRO',model.options.showprogress);

% FMINCON callbacks can be used, except that we must ensure the model is
% sent to the callbacks also (KNITRO only sends x)
funcs.objective = @(x)fmincon_fun(x,model);
funcs.constraints = @(x)fmincon_con(x,model);

switch model.options.verbose
    case 0
        model.options.knitro.Display = 'off';
    otherwise
        model.options.knitro.Display = 'iter';    
end

model.extendedFeatures = [];
model.objFnType = [];
model.xType = zeros(length(model.lb),1);
model.xType(model.binary_variables) = 2;
model.xType(model.integer_variables) = 1;
model.cineqFnType = repmat(2,length(model.bnonlinineq),1);

solvertime = clock;
[x,fval,exitflag,output,lambda] = knitromatlab_mip(funcs.objective,model.x0,model.A,model.b,model.Aeq,model.beq,model.lb,model.ub,funcs.constraints,model.xType,model.objFnType,model.cineqFnType,model.extendedFeatures,model.options.knitro);
solvertime = etime(clock,solvertime);

x = RecoverNonlinearSolverSolution(model,x);

% Internal format for duals
D_struc = [];

% Check, currently not exhaustive...
problem = 0;
switch exitflag
    case 0
        problem = 0;
    case {-200,-204,-515}
        problem = 1;    
    case -101
        problem = 2;
    case -400
        problem = 3;        
    otherwise
        problem = 11;
end
        
% Save all data sent to solver?
if model.options.savesolverinput
    solverinput.model = model;  
    solverinput.funcs = funcs;  
else
    solverinput = [];
end

% Save all data from the solver?
if model.options.savesolveroutput
    solveroutput.x = x;
    solveroutput.fmin = fmin;
    solveroutput.exitflag = exitflag;
    solveroutput.output=output;
    solveroutput.lambda=lambda;
else
    solveroutput = [];
end

% Standard interface
output = createoutput(x,D_struc,[],problem,'KNITRO',solverinput,solveroutput,solvertime);