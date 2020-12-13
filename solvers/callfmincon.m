function output = callfmincon(model)

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

switch model.options.verbose
    case 0
        model.options.fmincon.Display = 'off';   
    otherwise
        model.options.fmincon.Display = 'iter';
end

if isfield(model.options.fmincon,'LargeScale')
    if isequal(model.options.fmincon.LargeScale,'off')
        model.A = full(model.A);
        model.b = full(model.b);
        model.Aeq = full(model.Aeq);
        model.beq = full(model.beq);
    end
end
if isfield(model.options.fmincon,'Algorithm')
    if isequal(model.options.fmincon.Algorithm,'sqp') || isequal(model.options.fmincon.Algorithm,'active-set')
        model.A = full(model.A);
        model.b = full(model.b);
        model.Aeq = full(model.Aeq);
        model.beq = full(model.beq);
    end
end

if model.derivative_available
    model.options.fmincon.GradObj = 'on';  
    model.options.fmincon.GradConstr = 'on';
end

if model.options.savedebug
    ops = model.options.fmincon;
    save fmincondebug model %A b Aeq beq x0 lb ub ops
end

if strcmp(model.options.fmincon.Algorithm,'trust-region-reflective')
    if ~model.linearconstraints
        model.options.fmincon.Algorithm = 'interior-point';
    elseif ~isempty(model.A)
        model.options.fmincon.Algorithm = 'interior-point';
    elseif model.nonlinearinequalities>0 | model.nonlinearequalities>0
        model.options.fmincon.Algorithm = 'interior-point';
    elseif any(~isinf(model.lb) | ~isinf(model.ub)) & ~isempty(model.Aeq)
        model.options.fmincon.Algorithm = 'interior-point';
    end
end

% Figure out which variables are artificially introduced to normalize
% arguments in callback operators to simplify chain rules etc. We can do
% our function evaluations and gradient computations in our lifted world,
% but only expose the model in the original variables to the nonlinear
% solver. 
% model = compressLifted(model);

global latest_xevaled
global latest_x_xevaled
latest_xevaled = [];
latest_x_xevaled = [];

showprogress('Calling FMINCON',model.options.showprogress);

if model.linearconstraints
    g = [];
else
    g = @(x)fmincon_con_liftlayer(x,model);
end

solvertime = tic;
f = @(x)fmincon_fun_liftlayer(x,model);
if ~exist('OCTAVE_VERSION','builtin')
    [xout,fmin,flag,output,lambda] = fmincon(f,model.x0,model.A,model.b,model.Aeq,model.beq,model.lb,model.ub,g,model.options.fmincon);
else
    [xout,fmin,flag,output] = fmincon(f,model.x0,model.A,model.b,model.Aeq,model.beq,model.lb,model.ub,g,model.options.fmincon);
    lambda = [];
end
solvertime = toc(solvertime);

if ~isempty(xout) && ~isempty(model.lift);
    x = zeros(length(model.linearindicies),1);
    x(model.lift.linearIndex) = xout;
    x(model.lift.liftedIndex) = model.lift.T*xout + model.lift.d;
    x = RecoverNonlinearSolverSolution(model,x);
else
    x = RecoverNonlinearSolverSolution(model,xout);
end

problem = 0;

% Internal format for duals
D_struc = [];

% Check, currently not exhaustive...
if flag==0
    problem = 3;
else
    if flag>0
        problem = 0;
    else
        if flag == -2
            problem = 1;
        else
            if isempty(x)
                x = repmat(nan,length(model.c),1);
            end
            if model.c'*x<-1e10 % Likely unbounded
                problem = 2;
            else          % Probably convergence issues
                problem = 5;
            end
        end
    end
end

% Save all data sent to solver?
if model.options.savesolverinput
    solverinput.model = model;
    solverinput.A = model.A;
    solverinput.b = model.b;
    solverinput.Aeq = model.Aeq;
    solverinput.beq = model.beq;
    solverinput.options = model.options.fmincon;
else
    solverinput = [];
end

% Save all data from the solver?
if model.options.savesolveroutput
    solveroutput.x = x;
    solveroutput.fmin = fmin;
    solveroutput.flag = flag;
    solveroutput.output=output;
    solveroutput.lambda=lambda;
else
    solveroutput = [];
end

% Standard interface
output = createoutput(x,D_struc,[],problem,'FMINCON',solverinput,solveroutput,solvertime);