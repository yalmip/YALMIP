function sys = optimizer(Constraints,Objective,options,x,u)
%OPTIMIZER  Container for optimization problem
%
%   OPT = OPTIMIZER(Constraints,Objective,options,x,u) exports an object that
%   contains precompiled numerical data to be solved for varying arguments
%   x, returning the optimal value of the expression u.
%
%   OPTIMIZER works most efficiently if the varying data x enters the
%   optmization problem affinely. For the general cae, much more logic has
%   to be applied when instantiating the numerical data for a parametric
%   value, and when compiling the model, it is harder fpr YALMIP to
%   understand what kind of model it will be once th parameters are fixed.
%
%   By default, display is turned off (since optimizer is used in
%   situations where many problems are solved repeatedly. To turn on
%   display, set the verbose option in sdpsetting to 2.
%
%   Example
%
%    The following problem creates an LP with varying upper and lower
%    bounds on the decision variable.
%   
%     A = randn(10,3);
%     b = rand(10,1)*19;
%     c = randn(3,1);
%
%     z = sdpvar(3,1);
%     sdpvar UB LB
%
%     Constraints = [A*z <= b, LB <= z <= UB];
%     Objective = c'*z
%     % We want the optimal z as a function of [LB;UB]
%     optZ = optimizer(Constraints,Objective,[],[LB; UB],z);
%     
%     % Compute the optimal z when LB=1, UB = 3;
%     zopt = optZ([1; 3])
%
%     % Compute two solutions, one for (LB,UB) [1;3] and one for (LB,UB) [2;6]
%     zopt = optZ([[1; 3], [2;6]])
%
%     % A second output argument can be used to catch infeasibility
%     [zopt,infeasible] = optZ([1; 3])
%
%     % To avoid the need to vectorize in order to handle multiple
%     parameters, a cell-based definition can be used 
%
%     optZ = optimizer(Constraints,Objective,[],{LB,UB},{z,sum(z)})
%     [zopt,infeasible] = optZ({1,3});
%     zopt{1}
%     zopt{2}

if nargin < 5
    error('OPTIMIZER requires 5 inputs');
end

% With the new optional cell-based format, the internal format is always a
% vector with all information stacked, both in and out. Hence, we need to
% save original sizes before stacking things up
xoriginal = x;
if ~isa(x,'cell')
    xoriginal = {x};
end
if isa(x,'cell')
    xvec = [];
  for i = 1:length(x)
      if ~(isa(x{i},'sdpvar') | isa(x{i},'ndsdpvar'))
          error(['The parameters must be SDPVAR objects. Parameter #' num2str(i) ' is a ' upper(class(x{i}))]);
      end
      if is(x{i},'complex')
          x{i} = [real(x{i});imag(x{i})];
          complexInput(i) = 1;
      else
          complexInput(i) = 0;
      end
      sizeOrigIn{i} = size(x{i});
      z = x{i}(:);
      mask{i} = uniqueNonZeroRows(z);      
      xvec = [xvec;z(mask{i})];
  end
  x = xvec;
else   
    if ~isreal(x)%,'complex')
        complexInput(1) = 1;
        x = [real(x);imag(x)];
    else
        complexInput(1) = 0;
    end
    sizeOrigIn{1} = size(x);
    x = x(:);
    mask{1} = uniqueNonZeroRows(x);
    x = x(mask{1});
end
nIn = length(x);
mIn = 1;

if isa(x,'sdpvar')
    if ~is(x,'lpcone')
    error('All parameter arguments have to be simple variables (i.e., not expressions such a+b or 1+a)');
    end
end

if isa(u,'cell')
    uvec = []; 
    for i = 1:length(u)        
        if is(u{i},'complex')
            complexOutput(i) = 1;
            u{i} = [real(u{i});imag(u{i})];
        else
            complexOutput(i) = 0;
        end        
        sizeOrigOut{i} = size(u{i});
        uvec = [uvec;u{i}(:)];
    end
    u = uvec;
else
    if (isa(u,'sdpvar') || isa(u,'ndsdpvar')) && is(u,'complex')
        complexOutput(1) = 1;
        u = [real(u);imag(u)];
    else
        complexOutput(1) = 0;
    end
    sizeOrigOut{1} = size(u);
    u = u(:);
end
nOut = length(u);
mOut = 1;

if isempty(options)
    options = sdpsettings;
end

if ~isa(options,'struct')
    error('Third argument in OPTIMIZER should be an options structure.');
end

% Silent by default. If we want displays, set to 2
options.verbose = max(options.verbose-1,0);

% Since code is based on a fake equality, we must avoid bound propagation
% based on equalities
options.avoidequalitybounds=1;

% Normalize...
if isa(Constraints,'constraint')
    Constraints = lmi(Constraints);
end

if ~isempty(Constraints)
    if ~isa(Constraints,'constraint') &  ~isa(Constraints,'lmi')
        error('The first argument in OPTIMIZER should be a set of constraints');
    end
end

if ~isempty(Constraints)
    if any(is(Constraints,'sos'))
        tempOps = options;
        tempOps.sos.model = 2;
        tempOps.verbose = max(0,tempOps.verbose-1);
        parameter_sos = [x;u;recover(depends(Objective))];
        parameter_sos = depends(parameter_sos);
        for i = 1:length(Constraints)
            if ~is(Constraints,'sos')
                parameter_sos = [parameter_sos depends(Constraints(i))];
            end
        end
        parameter_sos = recover(parameter_sos);
        [Constraints,Objective] = compilesos(Constraints,Objective,tempOps,parameter_sos);
    end
end

if ~isequal(options.solver,'')
    % User has specified solver. Let us impose this solver forcefully to
    % the compilation code, in order to handle nonlinear parameterizations
    if ~strcmp(options.solver(1),'+')
        options.solver = ['+' options.solver];
    end
end

if options.removeequalities
    error('''removeequalities'' in optimizer objects not allowed.');
end

if ~isempty(Constraints) & any(is(Constraints,'uncertain'))
    [Constraints,Objective,failure] = robustify(Constraints,Objective,options);
    [aux1,aux2,aux3,model] = export((x == repmat(pi,nIn*mIn,1))+Constraints,Objective,options,[],[],0);
else
    [aux1,aux2,aux3,model] = export((x == repmat(pi,nIn*mIn,1))+Constraints,Objective,options,[],[],0);    
end

if ~isempty(aux3)
    if isstruct(aux3)
        if ismember(aux3.problem, [-9 -5 -4 -3 -2 -1 14])
            error(['Failed exporting the model: ' aux3.info])
        end
    end
end

if norm(model.F_struc(1:nIn*mIn,1)-repmat(pi,length(x),1),inf) > 1e-10
    error('Failed exporting the model (try to specify another solver)')        
end

% Try to set up an optimal way to compute the output
base = getbase(u);
if isempty(u) || (is(u,'linear') & all(sum(base | base,2) == 1) & all(sum(base,2)==1) & all(base(:,1)==0))
    % This is just a vector of variables
    z = [];
    map = [];
    uvec = u(:);
    % Setup to do fast getvariables(uvec(i))
    U = getbase(uvec);
    Uvar = getvariables(uvec);
    U = U(:,2:end);
    [ii,jj,ss] = find(U');
    for i = 1:length(uvec)
        %var = getvariables(uvec(i)); Slow
        %var = Uvar(find(U(i,:)));
        var = Uvar(ii(i));
        mapIndex = find(var == model.used_variables);
        if ~isempty(mapIndex)
            map = [map;mapIndex];
        else
            map = [map;0];
        end
    end
else
    % Some expression which we will use assign and double to evaluate
    vars = depends(u);
    z = recover(vars);    
    map = [];
    for i = 1:length(z)        
        var = vars(i);
        mapIndex = find(var == model.used_variables);
        if ~isempty(mapIndex)
            map = [map;mapIndex];
        else
            map = [map;0];
        end
    end        
end

if ~isempty(u) && (isempty(map) | min(size(map))==0)
    error('The requested decision variable (argument 4) is not in model');
end

model.getsolvertime = 0;

model.solver.callhandle = str2func(model.solver.call);

model.options = pruneOptions(model.options);
model.hashCache = gen_rand_hash(0,size(model.monomtable,2),1);

sys.recover = aux2;
sys.model = model;
sys.dimin = [nIn mIn];
sys.dimout = [nOut mOut];
sys.diminOrig = sizeOrigIn;
sys.dimoutOrig = sizeOrigOut;
sys.complexInput = complexInput;
sys.complexOutput = complexOutput;
sys.mask = mask;
sys.map = map;
sys.input.xoriginal = xoriginal;
sys.input.expression = x;
sys.output.expression = u;
sys.output.z = z;
sys.lastsolution = [];
sys.ParametricSolution = [];
sys.model.infeasible = 0;
% This is not guaranteed to give the index in the order the variables where
% given (tested in test_optimizer2
% [a,b,c] = find(sys.model.F_struc(1:prod(sys.dimin),2:end));
% Could be done using
[b,a,c] = find(sys.model.F_struc(1:prod(sys.dimin),2:end)');
% but let us be safe
%b = [];
%for i = 1:prod(sys.dimin)
%    b = [b;find(sys.model.F_struc(i,2:end))];
%end
sys.model.parameterIndex = b;
used_in = find(any(sys.model.monomtable(:,b),2));
Q = sys.model.Q;
Qa = Q;Qa(:,b)=[];Qa(b,:)=[];
Qb = Q(:,b);Qb(b,:)=[];
if nnz(Q)>0
    zeroRow = find(~any(Q,1));
    Qtest = Q;Q(zeroRow,:)=[];Q(:,zeroRow)=[];
    problematicQP = nonconvexQuadratic(Qtest);%min(eig(full(Qtest)))<-1e-14;
else
    problematicQP = 0;
end
if  any(sum(sys.model.monomtable(used_in,:),2)>1) | any(sum(sys.model.monomtable(used_in,:) | sys.model.monomtable(used_in,:),2) > 1) | problematicQP | ~isempty(sys.model.evalMap) | any(any(sys.model.monomtable<0))
    sys.nonlinear = 1;
else
    sys.nonlinear = 0;
end
sys.F = Constraints;
sys.h = Objective;
sys.ops = options;

sys.complicatedEvalMap = 0;
% Are all nonlinear operators acting on simple parameters? Elimination
% strategy will only be applied on simple problems such as x<=exp(par)
if ~model.solver.evaluation
    for i = 1:length(sys.model.evalMap)
        if ~all(ismember(sys.model.evalMap{i}.variableIndex,sys.model.parameterIndex))
            sys.complicatedEvalMap = 1;
        end
        if length(sys.model.evalMap{i}.arg)>2
            sys.complicatedEvalMap = 1;
        end
    end
    if sys.complicatedEvalMap
        error('Parameters are currently only allowed to enter function such as exp, sin etc as exp(a), sin(b) etc.')
    end
end

sys.model.evalParameters = [];
if sys.nonlinear
    % These artificial equalities are removed if we will use eliminate variables
 
    %   sys.model.F_struc(1:length(sys.parameters),:) = [];
    %   sys.model.K.f = sys.model.K.f - length(sys.parameters);
    
    % Which variables are simple nonlinear operators acting on parameters
    evalParameters = [];
    for i = 1:length(sys.model.evalMap)
        if all(ismember(sys.model.evalMap{i}.variableIndex,sys.model.parameterIndex))
            evalParameters = [evalParameters;sys.model.evalMap{i}.computes(:)];
        end
    end
    sys.model.evalParameters = evalParameters;
end

% In case we perform partial instantiation, we have to remember where we
% came from originallty when finally solving problems
sys.instatiatedvalues = zeros(length(model.used_variables),1);
sys.orginal_usedvariables = sys.model.used_variables;
sys.orginal_parameterIndex = sys.model.parameterIndex;

sys.input.stochastics = cell(1,length(sys.diminOrig));
if ~isempty(Constraints)
    randDefinitions = find(is(Constraints,'random'));
    if ~isempty(randDefinitions)
        
        for i = 1:length(randDefinitions)
            Fi = Constraints(randDefinitions(i));
            randDef{i}.distribution = struct(struct(Fi).clauses{1}.data).extra.distribution;
            randDef{i}.variables = sdpvar(Fi);
            
            for j = 1:length(sys.diminOrig)
                if isequal(getbase(sys.input.xoriginal{j}),getbase(randDef{i}.variables)) && isequal(getvariables(sys.input.xoriginal{j}),getvariables(randDef{i}.variables))
                    sys.input.stochastics{j} = randDef{i}.distribution;                
                end
            end
        end
    end
end

% Remove place holder constraints. No longer used
sys.model.F_struc(1:prod(sys.dimin),:)=[];
sys.model.K.f = sys.model.K.f-prod(sys.dimin);

sys = class(sys,'optimizer');
sys = optimizer_precalc(sys);

function i = uniqueNonZeroRows(x);
B = getbase(x);
% Quick check for trivially unique rows, typical 99% case
[n,m] = size(B);
if n == m-1 && nnz(B)==n
    if isequal(B,[spalloc(n,1,0) speye(n)])
        i = 1:n;
        return
    end
end
if  length(unique(B*randn(size(B,2),1))) == n && nnz(B*randn(size(B,2),1)) == size(B,1)
    i = 1:n;
    return
end
[temp,i] = unique(B*randn(size(B,2),1));
z = find(~any(B,2));
i = setdiff(i,z);
i = i(:);