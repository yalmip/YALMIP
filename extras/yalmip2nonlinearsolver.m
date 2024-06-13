function model = yalmip2nonlinearsolver(model)

global newmodel

newmodel = 1;

model.dense = 0;
if ~(model.equalitypresolved || model.presolved)
    model = propagate_bounds_from_equalities(model);
end

K = model.K;
lb = model.lb;
ub = model.ub;
x0 = model.x0;
c = model.c;

% Pick out the positive conditions from cones ||Ax+b|| <= c'*x+d which will
% be treated as (Ax+b)'*(ax+b) <= (c'*x+d)^2,  c'*x+d >= 0
model = bounds_from_cones_to_lp(model);

if isempty(model.evaluation_scheme)
    model = build_recursive_scheme(model);
end
model = compress_evaluation_scheme(model);

% Do some pre-calc to be used in calls from fmincon
nonlinearindicies = union(find(model.variabletype~=0),model.evalVariables);
linearindicies    = setdiff(find(model.variabletype==0),nonlinearindicies);
model.nonlinearindicies = nonlinearindicies;
model.linearindicies    = linearindicies;

model.Anonlinineq = [];
model.bnonlinineq = [];
model.Anonlineq = [];
model.bnonlineq = [];

% Extract linear and nonlinear equality constraints
if K.f>0
    Aeq = -model.F_struc(1:1:K.f,2:end);
    beq = model.F_struc(1:1:model.K.f,1);
    
    nonlinear_equalities_indicies = find(any(Aeq(:,nonlinearindicies),2));
    model.Anonlineq = Aeq(nonlinear_equalities_indicies,:);
    model.bnonlineq = beq(nonlinear_equalities_indicies);
    
    Aeq(nonlinear_equalities_indicies,:) = [];
    beq(nonlinear_equalities_indicies,:) = [];
    Aeq(:,nonlinearindicies) = [];
    model.F_struc(1:model.K.f,:) = [];
    model.K.f = 0;
else
    Aeq = [];
    beq = [];
end

% Find nonlinear eualities implied by lower and upper bounds
if ~isempty(ub) && ~isempty(lb)
    nonlinearequality = find(lb(nonlinearindicies) == ub(nonlinearindicies));
    if ~isempty(nonlinearequality)
        for i = 1:length(nonlinearequality)
          %  model.Anonlineq = [model.Anonlineq;eyev(length(c),nonlinearindicies(nonlinearequality(i)))'];
            model.Anonlineq = [model.Anonlineq;sparse(1,nonlinearindicies(nonlinearequality(i)),1,1,length(c))];
            model.bnonlineq = [model.bnonlineq;lb(nonlinearindicies(nonlinearequality(i)))];
        end
    end
end

% Extract linear and nonlinear inequality constraints
if any(model.K.l)
    A = -model.F_struc(1:model.K.l,2:end);
    b = model.F_struc(1:model.K.l,1);
    
    nonlinear_inequalities_indicies = find(any(A(:,nonlinearindicies),2));
    
    model.Anonlinineq = A(nonlinear_inequalities_indicies,:);
    model.bnonlinineq = b(nonlinear_inequalities_indicies);
    
    A(nonlinear_inequalities_indicies,:) = [];
    b(nonlinear_inequalities_indicies,:) = [];
    A(:,nonlinearindicies) = [];
    
    model.F_struc(1:model.K.l,:) = [];
    model.K.l = 0;
else
    A = [];
    b = [];
end

% This helps with robustness in bnb in some cases
x0candidate = zeros(length(c),1);
if ~isempty(lb) && ~isempty(ub)
    bounded = find(~isinf(lb) & ~isinf(ub));
    x0candidate(bounded) = (lb(bounded) + ub(bounded))/2;
    bounded_below = find(~isinf(lb) & isinf(ub));
    x0candidate(bounded_below) = lb(bounded_below) + 0.5;
    bounded_above = find(~isinf(lb) & isinf(ub));
    x0candidate(bounded_above) = lb(bounded_above) + 0.5;
end

if isempty(x0)
    x0 = x0candidate(linearindicies);
else
    if ~isempty(lb) && ~isempty(ub)
        x0((x0 < lb) | (x0 > ub)) = x0candidate((x0 < lb) | (x0 > ub));
    end
    x0 = x0(linearindicies);
end

if ~isempty(lb)
    lb = lb(linearindicies);
end
if ~isempty(ub)
    ub = ub(linearindicies);
end

lb_old = lb;
ub_old = ub;
[lb,ub,A,b] = remove_bounds_from_Ab(A,b,lb,ub);
[lb,ub,Aeq,beq] = remove_bounds_from_Aeqbeq(Aeq,beq,lb,ub);

if any(model.variabletype == 4)
    problematic = find(any(model.monomtable(:,linearindicies) < 0 ,1));
    if ~isempty(problematic)
        problematic = problematic(find(x0(problematic)==0));
        Oneisfeas = problematic(find(ub(problematic) > 1));
        x0(Oneisfeas) = 1;
    end
    
    problematic = find(any(model.monomtable(:,linearindicies)~=fix(model.monomtable(:,linearindicies)) ,1));
    lb(problematic) = max(lb(problematic),0);
end
x0(find(lb==ub)) = lb(find(lb==ub));
    
if size(A,1) == 0
    A = [];
end

if size(b,1) == 0
    b = [];
end

if size(Aeq,1) == 0
    Aeq = [];
end

if size(beq,1) == 0
    beq = [];
end

if model.presolveequalities
    if ~isempty(beq) &  (~model.equalitypresolved | ~(isequal(lb,lb_old) & isequal(ub,ub_old)))
        % This helps when there are artificial variables introduced to model
        % nonlinear operators such as log(2*x+1)
        p.F_struc = [beq -Aeq];
        p.K.f = size(beq,1);
        p.lb = lb;
        p.ub = ub;
        p.variabletype = zeros(1,length(lb));
        p.binary_variables = [];
        p.integer_variables = [];
        p = propagate_bounds_from_equalities(p);
        lb = p.lb;
        ub = p.ub;
    end
end

model.A = A;
model.b = b;
model.Aeq = Aeq;
model.beq = beq;
model.lb = lb;
model.ub = ub;
model.x0 = x0;

model = setup_fmincon_params(model);

% Check if all derivatives are available
model.derivative_available = 1;
for i = 1:length(model.evalMap)
    if isempty(model.evalMap{i}.properties.derivative)
        model.derivative_available = 0;
        break
    end
end

% Some precomputation of computational scheme for Jacobian
allA = [model.Anonlineq;model.Anonlinineq];
if anyCones(model.K)
    allA = [allA;model.F_struc(startofSOCPCone(model.K):end,2:end)];
end
requested = any(allA',2);
[i,j] = find((model.deppattern(find(requested),:)));
requested(j) = 1;
if ~isempty(model.evalMap)
    % Recursive stuff is only possible if we have evaluation-based
    % operators
    model.Crecursivederivativeprecompute = precomputeDerivative(model,requested);
end

% Some precomputation of computational scheme for gradient
requested = model.c | any(model.Q,2);
[i,j,k] = find((model.deppattern(find(requested),:)));
requested(j) = 1;
model.frecursivederivativeprecompute = precomputeDerivative(model,requested);

% Precomputed list of bilinear expressions, used in
% apply_recursive_differentiation
model = compile_bilinearslist(model);
model = compile_quadraticslist(model);

model.binary_variables  = find(ismember(linearindicies,model.binary_variables));
model.integer_variables  = find(ismember(linearindicies,model.integer_variables));
model.semicont_variables  = find(ismember(linearindicies,model.semicont_variables));
