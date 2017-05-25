function output = callfmincon(model)

model = yalmip2nonlinearsolver(model);

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
    callback_con = [];
else
    callback_con = 'fmincon_con_liftlayer';
end

solvertime = tic;
[xout,fmin,flag,output,lambda] = fmincon('fmincon_fun_liftlayer',model.x0,model.A,model.b,model.Aeq,model.beq,model.lb,model.ub,callback_con,model.options.fmincon,model);
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