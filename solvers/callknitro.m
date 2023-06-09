function output = callknitro(model)

% Extract and clear complementarity structure
if model.K.c(1) > 0
    C = model.F_struc(model.K.f+model.K.l+1:model.K.f+model.K.l+2*sum(model.K.c),:);
    model.F_struc(model.K.f+model.K.l+1:model.K.f+model.K.l+2*sum(model.K.c),:)=[];
    %model.K.l = model.K.l + sum(2*model.K.c);
else
    C = [];
end

% Standard NONLINEAR setup
tempF = model.F_struc;
tempK = model.K;
tempx0 = model.x0;
model = yalmip2nonlinearsolver(model);

%% Try to propagate initial values now that nonlinear evaluation logic is
%% set up. Do this if some of the intials are nan
model = propagateInitial(model,tempF,tempK,tempx0);

% Figure out which variables are artificially introduced to normalize
% arguments in callback operators to simplify chain rules etc. We can do
% our function evaluations and gradient computations in our lifted world,
% but only expose the model in the original variables to the nonlinear
% solver. 
if isempty(C) % Only do this if we don't have complementarity. FIXME
%    model = compressLifted(model);
end

if model.derivative_available
    model.options.knitro.GradObj = 'on';
    model.options.knitro.GradConstr = 'on';
else
    model.options.knitro.GradObj = 'off';
    model.options.knitro.GradConstr = 'off';
end
model.options.knitro.JacobPattern = jacobiansparsityfromnonlinear(model,0);

% If quadratic objective and no nonlinear constraints, we can supply an
% Hessian of the Lagrangian
usedinObjective = find(model.c | any(model.Q,2));
if ~any(model.variabletype(usedinObjective)) & any(any(model.Q))
    if ~any(model.K.q) && ~any(model.K.e) && ~any(model.K.p) && ~any(model.K.s)
        if  length(model.bnonlinineq)==0 & length(model.bnonlineq)==0
            H = model.Q(:,model.linearindicies);
            H = H(model.linearindicies,:);
            model.options.knitro.Hessian = 'user-supplied';
            model.options.knitro.HessPattern = sparse(H | H);
            model.options.knitro.HessFcn = @(x,l) 2*H;
        end
    end
end

global latest_xevaled
global latest_x_xevaled
global latest_x_g
global latest_x_f
global sdpLayer
latest_xevaled = [];
latest_x_xevaled = [];
latest_x_g = [];
latest_x_f = [];
sdpLayer.nullVectors = cell(length(model.K.s),1);
sdpLayer.eigenVectors = cell(length(model.K.s),1);
sdpLayer.oldGradient = cell(length(model.K.s),1);
sdpLayer.reordering  = cell(length(model.K.s),1);
sdpLayer.n  = inf;

if sdpLayer.n > max(model.K.s)
    sdpLayer.n  = inf;
end
if ~isinf(sdpLayer.n)
    % Prune jacobian sparsity
    J = model.options.knitro.JacobPattern(1:end-sum(model.K.s),:);
    Jsdp = model.options.knitro.JacobPattern(end-sum(model.K.s)+1:end,:);
    top = 1;
    for i = 1:length(model.K.s)
        J = [J;Jsdp(top:top + min(model.K.s(i),sdpLayer.n)-1,:)];
        top = top + model.K.s(i);
    end
    model.options.knitro.JacobPattern = J;
end

showprogress('Calling KNITRO',model.options.showprogress);

% FMINCON callbacks can be used, except that we 
% must ensure the model is sent to the callbacks also 
% KNITRO also calls twice, so we compute everything in first
% and then save stuff
funcs.constraints = @(x)knitro_callback_g(x,model);
funcs.objective = @(x)knitro_callback_f(x,model);

switch model.options.verbose
    case 0
        model.options.knitro.Display = 'off';
    otherwise
        model.options.knitro.Display = 'iter';
end

% SETUP complementarity information
if model.K.c(1) > 0
    top = 0;
    model.extendedFeatures.ccIndexList1 = [];
    model.extendedFeatures.ccIndexList2 = [];
    for i = 1:length(model.K.c)
        n = model.K.c(i);
        for j = 0:n-1
            j1 = find(C(top+i+j,:))-1;j1 = find(model.linearindicies == j1);
            j2 = find(C(top+i+n+j,:))-1;j2 = find(model.linearindicies == j2);
            model.lb(j1) = max(model.lb(j1),0);
            model.lb(j2) = max(model.lb(j2),0);
            model.extendedFeatures.ccIndexList1 = [model.extendedFeatures.ccIndexList1  j1];
            model.extendedFeatures.ccIndexList2 = [model.extendedFeatures.ccIndexList2 j2];
        end
        top = top + 2*n;
    end
else
    model.extendedFeatures = [];
end
model.objFnType = [];
model.xType = zeros(length(model.lb),1);
model.xType(model.binary_variables) = 2;
model.xType(model.integer_variables) = 1;
model.cineqFnType = repmat(2,length(model.bnonlinineq)+nnz(model.K.q),1);

solvertime = tic;
[xout,fval,exitflag,output,lambda] = knitromatlab_mip(funcs.objective,model.x0,model.A,full(model.b),model.Aeq,full(model.beq),model.lb,model.ub,funcs.constraints,model.xType,model.objFnType,model.cineqFnType,model.extendedFeatures,model.options.knitro,model.options.knitro.optionsfile);
solvertime = toc(solvertime);

if ~isempty(xout) && ~isempty(model.lift);
    x = zeros(length(model.linearindicies),1);
    x(model.lift.linearIndex) = xout(:);
    x(model.lift.liftedIndex) = model.lift.T*xout(:) + model.lift.d;
    x = RecoverNonlinearSolverSolution(model,x);
else
    x = RecoverNonlinearSolverSolution(model,xout);
end

% Internal format for duals
D_struc = [];

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
if model.options.savesolverinput
    solverinput.model = model;
    solverinput.funcs = funcs;
else
    solverinput = [];
end

% Save all data from the solver?
if model.options.savesolveroutput
    solveroutput.x = x;
    solveroutput.fval = fval;
    solveroutput.exitflag = exitflag;
    solveroutput.output=output;
    solveroutput.lambda=lambda;
else
    solveroutput = [];
end

% Standard interface
output = createOutputStructure(x,D_struc,[],problem,model.solver.tag,solverinput,solveroutput,solvertime);