function tests = test_model_review_regressions
tests = functiontests(localfunctions);

function test_scs_mixed_cones(testCase)
m = baseModel(1); m.K.e = 1; m.K.s = 2;
m.F_struc = sparse([[1;2;3;10;20;20;30] (1:7)']);
s = yalmip2scs(m);
testCase.verifyEqual(s.data.b,[10/sqrt(2);20;30/sqrt(2);1;2;3],'AbsTol',1e-12);
testCase.verifyEqual(full(s.data.A),-[4/sqrt(2);5;7/sqrt(2);1;2;3],'AbsTol',1e-12);

function test_copt_mixed_cones(testCase)
for linearRows = 0:2
    m = baseModel(2); m.K.q = 3; m.K.e = 1; m.K.l = linearRows;
    m.F_struc = sparse(linearRows+6,3);
    s = yalmip2copt(m);
    testCase.verifySize(s.A,[linearRows+6 8]);
    testCase.verifyEqual(full(s.A(:,3:5)),[zeros(linearRows,3);eye(3);zeros(3)]);
    testCase.verifyEqual(full(s.A(:,6:8)),[zeros(linearRows+3,3);flipud(eye(3))]);
end

function test_copt_bounds_only(testCase)
m = baseModel(2); m.lb = [0;0]; m.ub = [1;1];
s = yalmip2copt(m);
testCase.verifyEqual(size(s.A,1),numel(s.lhs));
testCase.verifyEqual(size(s.A,1),numel(s.rhs));
testCase.verifyEqual(size(s.A,2),2);

function test_gurobi_monomial_order(testCase)
for powers = {[1 2],[2 1],[1 3],[3 1],[2 3]}
    m = baseModel(3); m.c = [0;0;1];
    m.monomtable = sparse([1 0 0;0 1 0;powers{1} 0]);
    m.variabletype = [0 0 3]; m.solver.version = 'NONLINEAR';
    model = yalmip2gurobinonlinear(m);
    testCase.verifyGreaterThan(numel(model.obj),3);
    testCase.verifyEqual(size(model.Q,1),numel(model.obj));
end

function test_robust_norm_distinct_affine_terms(testCase)
x = sdpvar(1); w = sdpvar(1); Z.r = 1; Z.center = 0;
ops.verbose = 0; ops.robust.auxreduce = 'none'; V.w_variables = getvariables(w);
[F,~] = filter_normball([-w*(x+1)>=0,-w*(x-1)>=0],Z,x,w,w,'2-norm',ops,V);
vars = recover(depends(F)); assign(vars,zeros(length(vars),1)); assign(x,-1);
testCase.verifyLessThan(min(check(F)),0);

function test_polya_negative_coefficients(testCase)
w = sdpvar(2,1);
threw = false;
try
    filter_polya(lmi(-w(1)^2-w(2)^2>=0),w,0);
catch ME
    threw = ~isempty(strfind(ME.message,'Trivially infeasible'));
end
testCase.verifyTrue(threw);

function test_iff_parallel_inequalities(testCase)
x = sdpvar(1); b = binvar(1); bounds(x,-10,10);
cases = {[x<=1,x<=2],[x<=1,x>=2],[x<=2,x>=1]};
points = [3 1.5 3];
for i = 1:numel(cases)
    F = iff_internal(b,cases{i});
    aux = setdiff(depends(F),[getvariables(x),getvariables(b)]);
    assign(x,points(i)); assign(b,0);
    best = -Inf;
    for bits = 0:2^numel(aux)-1
        assign(recover(aux),bitget(bits,1:numel(aux))');
        best = max(best,min(check(F)));
    end
    testCase.verifyGreaterThanOrEqual(best,-1e-12);
end

function m = baseModel(n)
m.F_struc = sparse(0,n+1); m.c = zeros(n,1); m.Q = sparse(n,n); m.f = 0;
m.K.f = 0; m.K.l = 0; m.K.q = 0; m.K.e = 0; m.K.p = 0; m.K.s = 0;
m.K.sos.type = []; m.K.sos.variables = {}; m.K.sos.weight = {};
m.options.verbose = 0; m.options.scs = struct; m.options.copt = struct; m.options.gurobi = struct;
m.lb = -inf(n,1); m.ub = inf(n,1); m.x0 = [];
m.integer_variables = []; m.binary_variables = []; m.semicont_variables = []; m.evalMap = {};
