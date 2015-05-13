function output = callbonmin(model)

model = yalmip2nonlinearsolver(model);
options = [];
try
    options.bonmin = optiRemoveDefaults(model.options.bonmin,bonminset());
catch
    options.bonmin = model.options.bonmin;
end
options.ipopt = model.options.ipopt;
options.display = model.options.verbose;  

if ~model.derivative_available
    disp('Derivate-free call to bonmin/ipopt not yet implemented')
    error('Derivate-free call to bonmin/ipopt not yet implemented')
end

if model.options.savedebug
    save bonmindebug model
end

Fupp = [ repmat(0,length(model.bnonlinineq),1);
    repmat(0,length(model.bnonlineq),1);
    repmat(0,length(model.b),1);
    repmat(0,length(model.beq),1)];

Flow = [ repmat(-inf,length(model.bnonlinineq),1);
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
        output = createoutput(model.lb*0,[],[],problem,'BONMIN',solverinput,solveroutput,0);
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
latest_G= [];
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
    m = length(model.lb);    
    allA=[model.Anonlinineq; model.Anonlineq];
    jacobianstructure = spalloc(size(allA,1),m,0);    
    depends = allA | allA;   
    for i = 1:size(depends,1)
        vars = find(depends(i,:));
        [ii,vars] = find(model.deppattern(vars,:));
        vars = unique(vars);
        s = size(jacobianstructure,1);
        for j = 1:length(vars)            
            jacobianstructure(i,find(vars(j) == model.linearindicies)) = 1; 
        end      
    end
    allA=[model.A; model.Aeq];
    depends = allA | allA;
    jacobianstructure = [jacobianstructure;depends];
    
    Z = sparse(jacobianstructure);
    funcs.jacobianstructure = @() Z;
end

if ~model.options.usex0
    model.x0 = (options.lb+options.ub)/2;
    model.x0(isinf(options.ub)) = options.lb(isinf(options.ub))+1;
    model.x0(isinf(options.lb)) = options.ub(isinf(options.lb))-1;
    model.x0(isinf(model.x0)) = 0;
end

if ~isempty(model.binary_variables) | ~isempty(model.integer_variables)
    options.var_type = zeros(length(model.linearindicies),1);
    options.var_type(model.binary_variables) = -1;
    options.var_type(model.integer_variables) = 1;
end

showprogress('Calling BONMIN',model.options.showprogress);
solvertime = tic;
[xout,info] = bonmin(model.x0,funcs,options);
solvertime = toc(solvertime);

x = RecoverNonlinearSolverSolution(model,xout);

switch info.status
    case {0,2}
        problem = 0;
    case 1
        problem = 1;
    case {-1}
        problem = 3;
    case {3}
        problem = 15;
    otherwise
        problem = -1;
end

% Duals currently not supported
D_struc = [];

% Save all data sent to solver?
if model.options.savesolverinput
    solverinput.x0 = model.x0;
    solverinput.model = model;
    solverinput.options = options;
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
output = createoutput(x,D_struc,[],problem,'BONMIN',solverinput,solveroutput,solvertime);


% Code supplied by Jonatan Currie
function opts = removeDefaults(opts,defs)
oFn = fieldnames(opts);
for i = 1:length(oFn)
    label = oFn{i};
    if(isfield(defs,label))
        if(ischar(opts.(label)))
            if(strcmpi(defs.(label),opts.(label)))
                opts = rmfield(opts,label);
            end
        else
            if(defs.(label) == opts.(label))
                opts = rmfield(opts,label);
            end
        end
    end
end


