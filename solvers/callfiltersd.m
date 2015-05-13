function output = callfiltersd(model)

options = model.options.filtersd;

model = yalmip2nonlinearsolver(model);

% Make sure the callback returns a dense Jacobian
model.dense = 1;

if ~model.derivative_available
    disp('Derivate-free call to filterSD not yet implemented')
    error('Derivate-free call to filterSD not yet implemented')
end
if model.options.savedebug
    save filtersddebug model
end
showprogress('Calling filterSD',model.options.showprogress);

cu = [ repmat(0,length(model.bnonlinineq),1);
    repmat(0,length(model.bnonlineq),1);
    repmat(0,length(model.b),1);
    repmat(0,length(model.beq),1)];

cl = [ repmat(-inf,length(model.bnonlinineq),1);
    repmat(0,length(model.bnonlineq),1);
    repmat(-inf,length(model.b),1);
    repmat(0,length(model.beq),1)];

if isempty(cu)
    Flow = [];
    Fupp = [];
end

% These are needed to avoid recomputation due to ipopts double call to get
% f and df, and g and dg
global latest_x_f
global latest_x_g
global latest_df
global latest_f
global latest_G
global latest_g
global latest_xevaled
global latest_x_xevaled
latest_G = [];
latest_g = [];
latest_x_f = [];
latest_x_g = [];
latest_xevaled = [];
latest_x_xevaled = [];

funcs.objective = @(x)ipopt_callback_f(x,model);
funcs.gradient = @(x)ipopt_callback_df(x,model);
if ~isempty(cu)
    funcs.constraints = @(x)ipopt_callback_g(x,model);
    funcs.jacobian  = @(x)ipopt_callback_dg(x,model);
else
    funcs.constraints = [];
    funcs.jacobian = [];
end

lb  = model.lb(:)';
ub = model.ub(:)';

if ~model.options.usex0
    model.x0 = (lb+ub)/2;
    model.x0(isinf(ub)) = lb(isinf(ub))+1;
    model.x0(isinf(lb)) = ub(isinf(lb))-1;
    model.x0(isinf(model.x0)) = 0;
end

options.display = model.options.verbose;
solvertime = tic;
[xout,fval,exitflag,stats,lambda] = filtersd(funcs.objective, funcs.gradient, model.x0, lb, ub, funcs.constraints, funcs.jacobian, cl, cu, options);
solvertime = toc(solvertime);

% Duals currently not supported
lambda = [];

x = RecoverNonlinearSolverSolution(model,xout);

switch exitflag
    case {0}
        problem = 0;
    case {1}
        problem = 2;
    case {2,3}
        problem = 1;
    case {5}
        problem = 3;
    case {4,9}
        problem = 4;
    case 105
        problem = 16;        
    otherwise
        problem = -1;
end

% Internal format for duals
D_struc = [];

% Save all data sent to solver?
if model.options.savesolverinput
    solverinput.model = model;
else
    solverinput = [];
end

% Save all data from the solver?
if model.options.savesolveroutput
    solveroutput.x = xout;  
    solveroutput.exitflag = exitflag;
    solveroutput.stats = stats;
else
    solveroutput = [];
end

% Standard interface
output = createoutput(x,D_struc,[],problem,'filterSD',solverinput,solveroutput,solvertime);


