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
%   parameter (see Wiki for beta support of a much more general optimizer)
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
%   display, set the verbose option in sdpsetting to 2.
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
%     Constraints = [A*z <= b, LB <= z <= UB];
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
%     % A second output argument can be used to catch infeasibility
%     [zopt,infeasible] = optZ{[1; 3]}
%
%     % To avoid the need to vectorize in order to handle multiple
%       parameters, a cell-based format can be used, both for inputs and
%       outputs. Note that the optimizer object now is called with a cell
%       and returns a cell
%
%     optZ = optimizer(Constraints,Objective,[],{LB,UB},{z,sum(z)})
%     [zopt,infeasible] = optZ{{1,3}};
%     zopt{1}
%     zopt{2}

if nargin < 5
    error('OPTIMIZER requires 5 inputs');
end

% With the new optional cell-based format, the internal format is always a
% vector with all information stacked, both in and out. Hence, we need to
% save original sizes before stacking things up

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
      mask{i} = uniqueRows(z);      
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
    mask{1} = uniqueRows(x);
    x = x(mask{1});
end
nIn = length(x);
mIn = 1;

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
    if is(u,'complex')
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
if is(u,'linear') & all(sum(base | base,2) == 1) & all(sum(base,2)==1) & all(base(:,1)==0)
    % This is just a vecotr of variables
    z = [];
    map = [];
    uvec = u(:);
    for i = 1:length(uvec)
        var = getvariables(uvec(i));
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

if isempty(map) | min(size(map))==0
    error('The requested decision variable (argument 4) is not in model');
end

model.getsolvertime = 0;

model.solver.callhandle = str2func(model.solver.call);

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
sys.input.expression = x;
sys.output.expression = u;
sys.output.z = z;
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
sys.parameters = b;
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
for i = 1:length(sys.model.evalMap)
    if ~all(ismember(sys.model.evalMap{i}.variableIndex,sys.parameters))
       sys.complicatedEvalMap = 1;
    end
    if length(sys.model.evalMap{i}.arg)>2
        sys.complicatedEvalMap = 1;
    end
end

if sys.nonlinear & ~sys.complicatedEvalMap
    % These artificial equalities are removed if we will use eliminate variables
    sys.model.F_struc(1:length(sys.parameters),:) = [];
    sys.model.K.f = sys.model.K.f - length(sys.parameters);
    
    % Which variables are simple nonlinear operators acting on parameters
    evalParameters = [];
    for i = 1:length(sys.model.evalMap)
        if all(ismember(sys.model.evalMap{i}.variableIndex,sys.parameters))
            evalParameters = [evalParameters;sys.model.evalMap{i}.computes(:)];
        end
    end
    sys.model.evalParameters = evalParameters;
end

% This data is used in eliminatevariables (nonlinear parameterizations)
% A lot of performance is gained by precomputing them
% This will work as long as there a no zeros in the parameters, which might
% cause variables to dissapear (as in x*parameter >=0, parameter = 0)
% (or similiar effects happen)
sys.model.precalc.newmonomtable = sys.model.monomtable;
sys.model.precalc.rmvmonoms = sys.model.precalc.newmonomtable(:,sys.parameters);
sys.model.precalc.newmonomtable(:,sys.parameters) = 0;
sys.model.precalc.Qmap = [];
% R2012b...
try
    [ii,jj,kk] = stableunique(sys.model.precalc.newmonomtable*gen_rand_hash(0,size(sys.model.precalc.newmonomtable,2),1));
    sys.model.precalc.S = sparse(kk,1:length(kk),1);
    sys.model.precalc.skipped = setdiff(1:length(kk),jj);    
    sys.model.precalc.blkOneS = blkdiag(1,sys.model.precalc.S');     
catch  
end
    
if sys.nonlinear & ~sys.complicatedEvalMap
    
    % Precompute some structures
    newmonomtable = sys.model.monomtable;
    rmvmonoms = newmonomtable(:,[sys.parameters;evalParameters]);
    % Linear indexation to fixed monomial terms which have to be computed
    % [ii1,jj1] = find((rmvmonoms ~= 0) & (rmvmonoms ~= 1));
    [ii1,jj1] = find( rmvmonoms < 0 | rmvmonoms > 1 | fix(rmvmonoms) ~= rmvmonoms);    
    sys.model.precalc.index1 = sub2ind(size(rmvmonoms),ii1,jj1);    
    sys.model.precalc.jj1 = jj1;    
            
    % Linear indexation to linear terms
    linterms = rmvmonoms == 1;
    if ~isempty(jj1) | any(sum(linterms,2)>1)
        [ii2,jj2] = find(linterms);
        sys.model.precalc.index2 = sub2ind(size(rmvmonoms),ii2,jj2);
        sys.model.precalc.jj2 = jj2;
        sys.model.precalc.aux = rmvmonoms*0+1;
    else
        [ii2,jj2] = find(linterms);
        sys.model.precalc.index2 = ii2;
        sys.model.precalc.jj2 = jj2;
        sys.model.precalc.aux = ones(size(rmvmonoms,1),1);
    end
    
    sys.model.newmonomtable = model.monomtable;
    sys.model.rmvmonoms =  sys.model.newmonomtable(:,[sys.parameters;evalParameters]);
    sys.model.newmonomtable(:,union(sys.parameters,evalParameters)) = 0;
   
    sys.model.removethese = find(~any(sys.model.newmonomtable,2));
    sys.model.keepingthese = find(any(sys.model.newmonomtable,2));    
end

sys = class(sys,'optimizer');

function i = uniqueRows(x);
B = getbase(x);
% Quick check for trivially unique rows, typical 99% case
[n,m] = size(B);
if n == m-1 && nnz(B)==n
    if isequal(B,[spalloc(n,1,0) speye(n)])
        i = 1:n;
        return
    end
end
if  length(unique(B*randn(size(B,2),1))) == n
    i = 1:n;
    return
end
[temp,i,j] = unique(B,'rows');
i = i(:);