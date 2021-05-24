function output = callipoptgp(interfacedata)

% Retrieve needed data
options = interfacedata.options;
F_struc = interfacedata.F_struc;
c       = interfacedata.c;
Q       = interfacedata.Q;
K       = interfacedata.K;
x0      = interfacedata.x0;
extended_variables = interfacedata.extended_variables;
ub      = interfacedata.ub;
lb      = interfacedata.lb;
mt      = interfacedata.monomtable;

options = [];
try
    options.ipopt = optiRemoveDefaults(interfacedata.options.ipopt,ipoptset());
catch
    options.ipopt = options.ipopt;
end
options.ipopt.print_level = 2+interfacedata.options.verbose;

% *********************************
% What type of variables do we have
% *********************************
linear_variables = find((sum(abs(mt),2)==1) & (any(mt==1,2)));
nonlinear_variables = setdiff((1:size(mt,1))',linear_variables);
sigmonial_variables = find(any(0>mt,2) | any(mt-fix(mt),2));

% Convert to common format for fmincon, mosek and gpposy
ubtemp = ub;
lbtemp = lb;
fixed = find(ub(linear_variables) == lb(linear_variables));
ubtemp(linear_variables(fixed)) = inf;
lbtemp(linear_variables(fixed)) = -inf;
[prob,problem] = yalmip2geometric(options,F_struc,c,Q,K,ubtemp,lbtemp,mt,linear_variables,extended_variables);

%Add equalities for fixed variables
fixed_in_bound = [];
fixed_bound=[];
if ~isempty(fixed)
    prob.G = [prob.G;-sparse(1:length(fixed),fixed,1,length(fixed),length(linear_variables))];
    prob.h = [prob.h;lb(linear_variables(fixed))];
end

%something failed
if problem
    % Go to standard fmincon
    if interfacedata.options.verbose
        disp('Conversion to geometric program failed. Trying general non-convex model in ipopt');
        disp(' ');
    end
    interfacedata.solver.tag = strrep(interfacedata.solver.tag,'-geometric','');
    output = callipopt(interfacedata);
    return
end

if isempty(x0)
    x0 = zeros(length(linear_variables),1);
else
    x0 = x0(linear_variables);
end

% Fake logarithm (extend linearly for small values)
ind = find(x0<1e-2);
x0(ind) = exp(log(1e-2)+(x0(ind)-1e-2)/1e-2);
x0 = log(x0);

% Clean up the bounds (from branch and bound)
% Note, these bounds are in the
% logarithmic variables.
if ~isempty(lb)
    lb = lb(linear_variables);
    ind = find(lb<1e-2);
    lb(ind) = exp(log(1e-2)+(lb(ind)-1e-2)/1e-2);
    lb = log(lb+sqrt(eps));
end
if ~isempty(ub)
    ub = ub(linear_variables);
    ind = find(ub<1e-2);
    ub(ind) = exp(log(1e-2)+(ub(ind)-1e-2)/1e-2);
    ub = log(ub+sqrt(eps));
end

if interfacedata.options.savedebug
    ops = options.ipopt;
    save ipoptgpdebug prob x0 ops lb ub
end

prob = precalcgpstruct(prob);

Fupp = [repmat(0,size(prob.G,1)+max(prob.map),1);];
Flow = [repmat(-inf,max(prob.map),1);repmat(0,size(prob.G,1),1)];

% These are needed to avoid recomputation due to ipopts double call to get
% f and df, and g and dg (not fully done yet in geometric version)
global latest_x_f
global latest_x_g
global latest_df
global latest_f
global latest_G
global latest_g
latest_G = [];
latest_g = [];
latest_x_f = [];
latest_x_g = [];

funcs.objective = @(x)ipoptgp_callback_f(x,prob);
funcs.gradient = @(x)ipoptgp_callback_df(x,prob);
if ~isempty(Fupp)
    funcs.constraints = @(x)ipoptgp_callback_g(x,prob);
    funcs.jacobian  = @(x)ipoptgp_callback_dg(x,prob);
end

options.lb = lb(:)';
options.ub = ub(:)';
if ~isempty(Fupp)
    options.cl = Flow;
    options.cu = Fupp;
end

if ~isempty(Fupp)
    Z = double([(prob.B~=0)*(prob.A~=0);prob.G] | [prob.B*prob.A;prob.G]);
    Z = sparse(Z);
    funcs.jacobianstructure = @() Z;
end

solvertime = tic;
[xout,info] = ipopt(x0,funcs,options);
solvertime = toc(solvertime);

x = zeros(length(c),1);
x(linear_variables) = exp(xout);
x(fixed_in_bound) = fixed_bound;

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

% Internal format for duals, not supported yet
D_struc = [];

if interfacedata.options.savesolverinput
    solverinput.A = [];
    solverinput.b = [];
    solverinput.Aeq = [];
    solverinput.beq = [];
    solverinput.options = options;
else
    solverinput = [];
end

% Save all data from the solver?
if interfacedata.options.savesolveroutput
    solveroutput.xout = xout;    
    solveroutput.info = info;
else
    solveroutput = [];
end

% Standard interface
output = createOutputStructure(x,D_struc,[],problem,interfacedata.solver.tag,solverinput,solveroutput,solvertime);
