function output = callnomad(model)

%model.presolveequalities = 1;
%model.equalitypresolved = 0;
model = yalmip2nonlinearsolver(model);

% Nomad does not need derivatives, so let us inform our callbacks that we
% don't need there
model.derivative_available = 0;

% Note, nomad does not support equalities, so we place in equalities
nlrhs = [model.bnonlinineq*0;model.b*0;
         model.bnonlineq*0;model.bnonlineq*0;
         model.beq*0;model.beq*0];
if isempty(nlrhs)
    % clean
    nlrhs = [];
end
model.Anonlinineq = [model.Anonlinineq;model.Anonlineq;-model.Anonlineq];
model.Anonlineq = [];
model.bnonlinineq = [model.bnonlinineq;model.bnonlineq;-model.bnonlineq];
model.bnonlineq = [];
model.A = [model.A;-model.Aeq;model.Aeq];
model.Aeq = [];
model.b = [model.b;-model.beq;model.beq];
model.beq = [];

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
funcs.constraints = @(x)ipopt_callback_g(x,model);

lb = model.lb(:);
ub = model.ub(:);

if ~model.options.usex0
    model.x0 = (lb+ub)/2;
    model.x0(isinf(ub)) = lb(isinf(ub))+1;
    model.x0(isinf(lb)) = ub(isinf(lb))-1;
    model.x0(isinf(model.x0)) = 0;
end

integer_variables = model.integer_variables;
binary_variables = model.binary_variables;
xtype = char(ones(length(model.c),1)*67);
xtype(integer_variables) = 'I';
xtype(binary_variables)  = 'B';  % Should not happen except from bmibnb
xtype = xtype( model.linearindicies);
opts = model.options.nomad;
opts.display_degree = model.options.verbose;

showprogress('Calling NOMAD',model.options.showprogress);
solvertime = tic;
[x,fval,exitflag,iter,nfval] = nomad(funcs.objective,model.x0,lb,ub,funcs.constraints,nlrhs,xtype,opts);
solvertime = toc(solvertime);

% Duals currently not supported
lambda = [];

x = RecoverNonlinearSolverSolution(model,x);

switch exitflag
    case {1}
        problem = 0;
    case {-1}
        problem = 1;
    case {0}
        problem = 3;
    case {-2}        
        problem = 11;
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
    solveroutput.lambda = lambda;
    solveroutput.iters = iters;
    solveroutput.info = info;
else
    solveroutput = [];
end

% Standard interface
output = createoutput(x,D_struc,[],problem,'NOMAD',solverinput,solveroutput,solvertime);


