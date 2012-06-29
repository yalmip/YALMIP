function sys = optimizer(Constraints,Objective,options,x,u)
%OPTIMIZER  Container for optimization problem
%
%   OPT = OPTIMIZER(Constraints,Objective,options,x,u) exports an object that
%   contains precompiled numerical data to be solved for varying arguments
%   x, returning the optimal value of the expression u.
%
%   OPTIMIZER typically only works efficiently if the varying data x enters
%   the optmization problem affinely. If not, the precompiled problems will
%   be nonconvex, despite the problem being simple for a fixed value of the
%   parameter.
%
%   In principle, if an optimization problem has a parameter x, and we
%   repeatedly want to solve the problem for varying x to compute a
%   variable u, we can, instead of repeatedly constructing optimization
%   problems for fixed values of x, introduce a symbolic x, and then
%   simply add an equality 
%    solvesdp([Constraints,x == value],Objective);
%    uopt = double(u)
%   There will still be overhead from the SOLVESDP call, so we can
%   precompile the whole structure, and let YALMIP handle the addition of
%   the equality constraint for the fixed value, and automatically extract
%   the solution variables we are interested in
%    OPT = optimizer(Constraints,Objective,options,x,u)
%    uopt1 = OPT{value1}
%    uopt2 = OPT{value2}
%    uopt3 = OPT{value3}
%
%   By default, display is turned off (since optimizer is used in
%   situations where many problems are solved repeatedly. To turn on
%   display, set the vcerbose option in sdpsetting to 2.
%
%   Example
%
%    The following problem creates an LP with varying upper and lower
%    bounds on the decision variable.
%
%    The optimizing argument is obtained by indexing (with {}) the optimizer 
%    object with the point of interest. The argument should be a column
%    vector (if the argument has a width larger than 1, YALMIP assumes that
%    the optimal solution should be computed in several points) 
%   
%     A = randn(10,3);
%     b = rand(10,1)*19;
%     c = randn(3,1);
%
%     z = sdpvar(3,1);
%     sdpvar UB LB
%
%     Constraints = [A*z < b, LB < z < UB];
%     Objective = c'*z
%     % We want the optimal z as a function of [LB;UB]
%     optZ = optimizer(Constraints,Objective,[],[LB; UB],z);
%     
%     % Compute the optimal z when LB=1, UB = 3;
%     zopt = optZ{[1; 3]}
%
%     % Compute two solutions, one for (LB,UB) [1;3] and one for (LB,UB) [2;6]
%     zopt = optZ{[[1; 3], [2;6]]}
%
%     A second output argument can be used to catch infeasibility
%     [zopt,infeasible] = optZ{[1; 3]}

if nargin < 5
    error('OPTIMIZER requires 5 inputs');
end

[n,m] = size(x);
x = x(:);

if isempty(options)
    options = sdpsettings;
end
% Silent by default. If we want displays, set to 2
options.verbose = max(options.verbose-1,0);

% Since code is based on a fake equality, we must avoid bound propagation
% based on equalities
options.avoidequalitybounds=1;

% Normalize...
if isa(Constraints,'constraint')
    Constraints = set(Constraints);
end

if any(is(Constraints,'uncertain'))
    [Constraints,Objective,failure] = robustify(Constraints,Objective,options);
    [aux1,aux2,aux3,model] = export(set(x == repmat(pi,n*m,1))+Constraints,Objective,options,[],[],0);
else
    [aux1,aux2,aux3,model] = export(set(x == repmat(pi,n*m,1))+Constraints,Objective,options,[],[],0);    
end

if ~isempty(aux3)
    if isstruct(aux3)
        if aux3.problem == -4
            error(['Failed exporting the model: ' aux3.info])
        end
    end
end

if norm(model.F_struc(1:n,1)-repmat(pi,length(x),1),inf) > 1e-10
    error('Failed exporting the model (try to specify another solver)')    
end

% Try to set up an optimal way to compute the output
base = getbase(u);
if is(u,'linear') & all(sum(base | base,2) == 1) & all(sum(base,2))==1 & all(base(:,1)==0)
    % This is just a vecotr of variables
    z = [];
    map = [];
    for i = 1:length(u)
        var = getvariables(u(i));
        map = [map;find(var == model.used_variables)];
    end
else
    % Some expression which we will use assign and double to evaluate
    z = recover(depends(u));
    map = [];
    for i = 1:length(z)
        var = getvariables(z(i));
        map = [map;find(var == model.used_variables)];
    end        
end

if isempty(map) | min(size(map))==0
    error('The requested decision variable (argument 4) is not in model');
end

model.getsolvertime = 0;

model.solver.callhandle = str2func(model.solver.call);

sys.recover = aux2;
sys.model = model;
sys.dimin = [n m];
sys.dimout = size(u);
sys.map = map;
sys.input.expression = x;
sys.output.expression = u;
sys.output.z = z;
[a,b,c] = find(sys.model.F_struc(1:prod(sys.dimin),2:end));
sys.parameters = b;
used_in = find(any(sys.model.monomtable(:,b),2));
if any(sum(sys.model.monomtable(used_in,:) | sys.model.monomtable(used_in,:),2) > 1)
    sys.nonlinear = 1;
else
    sys.nonlinear = 0;
end
sys.F = Constraints;
sys.h = Objective;
sys.ops = options;
sys = class(sys,'optimizer');

