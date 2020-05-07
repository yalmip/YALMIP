function output = callipopt(model)

options = [];
try
    options.ipopt = optiRemoveDefaults(model.options.ipopt,ipoptset());
catch
    options.ipopt = model.options.ipopt;
end
options.ipopt.print_level = 2+model.options.verbose;

% Standard NONLINEAR setup
tempF = model.F_struc;
tempK = model.K;
tempx0 = model.x0;
model = yalmip2nonlinearsolver(model);
% % Try to propagate initial values now that nonlinear evaluation logic is
% % set up. Do this if some of the intials are nan
if any(isnan(model.x0)) & ~all(isnan(model.x0))        
    [~,~,xevaledout] = fmincon_fun(model.x0,model);
    startNan = nnz(isnan(xevaledout));
    goon = 1;
    while goon
        temp = propagatex0(xevaledout,tempF,tempK);
        model.x0 = temp(model.linearindicies);
        [~,~,xevaledout] = fmincon_fun(model.x0,model);
        goon =  nnz(isnan(xevaledout)) < startNan;
        startNan = nnz(isnan(xevaledout));
    end     
end
model.x0(isnan(model.x0))=0;


if ~model.derivative_available
    disp('Derivate-free call to ipopt not yet implemented')
    error('Derivate-free call to ipopt not yet implemented')
end
if model.options.savedebug
    save ipoptdebug model
end
showprogress('Calling IPOPT',model.options.showprogress);

% Figure out which variables are artificially introduced to normalize
% arguments in callback operators to simplify chain rules etc. We can do
% our function evaluations and gradient computations in our lifted world,
% but only expose the model in the original variables to the nonlinear
% solver. 
% model = compressLifted(model);

Fupp = [ repmat(0,length(model.bnonlinineq)+length(model.K.q)*(model.K.q(1)>0),1);
    repmat(0,length(model.bnonlineq),1);
    repmat(0,length(model.b),1);
    repmat(0,length(model.beq),1)];

Flow = [ repmat(-inf,length(model.bnonlinineq)+length(model.K.q)*(model.K.q(1)>0),1);
    repmat(0,length(model.bnonlineq),1);
    repmat(-inf,length(model.b),1);
    repmat(0,length(model.beq),1)];

if isempty(Flow)
    Flow = [];
    Fupp = [];
end

% Since ipopt react strangely on lb>ub, we should bail if that is detected
% (ipopt creates an exception)
if ~isempty(model.lb)
    if any(model.lb>model.ub)
        problem = 1;   
        solverinput = [];
        solveroutput = [];  
        output = createoutput(model.c*0,[],[],problem,'IPOPT',solverinput,solveroutput,0);
        return
    end
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
if ~isempty(Fupp)
    funcs.constraints = @(x)ipopt_callback_g(x,model);
    funcs.jacobian  = @(x)ipopt_callback_dg(x,model);
end

options.lb = model.lb(:)';
options.ub = model.ub(:)';
if ~isempty(Fupp)
    options.cl = Flow;
    options.cu = Fupp;
end

if ~isempty(Fupp)
    Z = jacobiansparsityfromnonlinear(model);
    funcs.jacobianstructure = @() Z;
end

if ~model.options.usex0
    model.x0 = (options.lb+options.ub)/2;
    model.x0(isinf(options.ub)) = options.lb(isinf(options.ub))+1;
    model.x0(isinf(options.lb)) = options.ub(isinf(options.lb))-1;
    model.x0(isinf(model.x0)) = 0;
    if any(model.variabletype == 4)
        problematic = find(any(model.monomtable(:,model.linearindicies) < 0 ,1));
        if ~isempty(problematic)
            problematic = problematic(find(model.x0(problematic)==0));
            Oneisfeas = problematic(find(model.ub(problematic) > 1));
            model.x0(Oneisfeas) = 1;
        end
    end
    model.x0(find(model.lb==model.ub)) = model.lb(find(model.lb==model.ub));
end

% If quadratic objective and no nonlinear constraints, we can supply an
% Hessian of the Lagrangian
usedinObjective = find(model.c | any(model.Q,2));
if ~any(model.variabletype(usedinObjective)) & any(model.Q)
    if  length(model.bnonlinineq)==0 & length(model.bnonlineq)==0
        H = model.Q(:,model.linearindicies);
        H = H(model.linearindicies,:);
        funcs.hessian = @(x,s,l) tril(2*H);
        funcs.hessianstructure = @()tril(sparse(double(H | H)));      
        options.ipopt.hessian_approximation = 'exact';
    end
end

solvertime = tic;
[xout,info] = ipopt(model.x0,funcs,options);
solvertime = toc(solvertime);

% Duals currently not supported
lambda = [];

if ~isempty(xout) && ~isempty(model.lift);
    x = zeros(length(model.linearindicies),1);
    x(model.lift.linearIndex) = xout(:);
    x(model.lift.liftedIndex) = model.lift.T*xout(:) + model.lift.d;
    x = RecoverNonlinearSolverSolution(model,x);
else
    x = RecoverNonlinearSolverSolution(model,xout);
end

switch info.status
    case {0,1}
        problem = 0;
    case {2}
        problem = 1;
    case {-1}
        problem = 3;
    case {3,4,-2,-3}
        problem = 4;
    case {-11,-12,-13}
        problem = 7;
    case {-10,-100,-101,-102,-199}
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
    solveroutput.info = info;
else
    solveroutput = [];
end

% Standard interface
output = createoutput(x,D_struc,[],problem,'IPOPT',solverinput,solveroutput,solvertime);


