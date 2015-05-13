function output = callopenopt(model)

% Retrieve needed data
options = model.options;
F_struc = model.F_struc;
c       = model.c;
K       = model.K;
x0      = model.x0;
Q       = model.Q;
lb      = model.lb;
ub      = model.ub;
monomtable = model.monomtable;

if isempty(model.evaluation_scheme)
    model = build_recursive_scheme(model);
end

nonlinearindicies = union(find(model.variabletype~=0),model.evalVariables);
linearindicies    = setdiff(find(model.variabletype==0),nonlinearindicies);
model.linearindicies    = linearindicies;
[lb,ub] = findulb(model.F_struc,model.K);

if isempty(model.x0)
    model.x0 = zeros(length(linearindicies),1);
end

model.Anonlinineq = [];
model.bnonlinineq = [];
model.Anonlineq = [];
model.bnonlineq = [];

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

if ~isempty(ub) & ~isempty(lb)
    nonlinearequality = find(lb(nonlinearindicies) == ub(nonlinearindicies));
    if ~isempty(nonlinearequality)
        for i = 1:length(nonlinearequality)
            model.Anonlineq = [model.Anonlineq;eyev(length(c),nonlinearindicies(nonlinearequality(i)))'];
            model.bnonlineq = [model.bnonlineq;lb(nonlinearindicies(nonlinearequality(i)))];
        end
    end
end

if model.K.l>0
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

if isfield(options.fmincon,'LargeScale')
    if isequal(options.fmincon.LargeScale,'off')
        A = full(A);
        b = full(b);
        Aeq = full(Aeq);
        beq = full(beq);
    end
end

lb = lb(linearindicies);
ub = ub(linearindicies);

% Initialize persistent variables for callback
openopt_fun([],model);

prob = ooAssign(@openopt_fun, model.x0,'MaxIter=1500');
prob.iterPrint = options.verbose-1;
prob.doPlot = 0;
prob.A = A;
prob.b = b;
prob.Aeq = Aeq;
prob.beq = beq;
prob.lb = lb;
prob.ub = ub;

solvertime = tic;
r = ooRun(prob, 'ralg');
solvertime = toc(solvertime);

if isempty(nonlinearindicies)
    x = r.xf(:);
else
    x = zeros(length(c),1);
    for i = 1:length(linearindicies)
        x(linearindicies(i)) = r.xf(i);
    end
    x = x(1:length(c));
end

switch r.istop
    case 4
        problem = 0;
    otherwise
            problem = 1;
end

% Internal format for duals
D_struc = [];

% Save all data sent to solver?
if options.savesolverinput
    solverinput = [];
else
    solverinput = [];
end

% Save all data from the solver?
if options.savesolveroutput
    solveroutput = [];
else
    solveroutput = [];
end

% Standard interface
output = createoutput(x,D_struc,[],problem,'SIMANN',solverinput,solveroutput,solvertime);