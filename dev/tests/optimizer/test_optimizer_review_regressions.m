function tests = test_optimizer_review_regressions
tests = functiontests(localfunctions);

function setupOnce(testCase)
testCase.TestData.options = sdpsettings('solver','+none','verbose',0);

function test_merge_equalities(testCase)
p = sdpvar(1); x = sdpvar(1); y = sdpvar(1);
ops = testCase.TestData.options;
a = optimizer([x==p,y>=0],x+y,ops,p,[x;y]);
b = optimizer([x==p,y==1,y>=0],x+y,ops,p,[x;y]);
c = [a b]; m = c.model;
testCase.verifyEqual(size(m.F_struc,1),m.K.f+m.K.l+sum(m.K.q)+sum(m.K.s.^2));
mb = b.model;
testCase.verifyEqual(m.K.f,mb.K.f);

function test_merge_sdp_blocks(testCase)
p = sdpvar(1); x = sdpvar(1); y = sdpvar(1);
ops = testCase.TestData.options;
F1 = [x 1;1 y]>=0; F2 = [x-5 2;2 y-4]>=0; F3 = [x-5 2;2 y-100]>=0;
a = optimizer([F1,x+y<=p],x+y,ops,p,[x;y]);
b = optimizer([F1,F2,x+y<=p],x+y,ops,p,[x;y]);
c = [a b];
mc = c.model; mb = b.model;
testCase.verifyEqual(mc.K.s,mb.K.s);
testCase.verifyEqual(mc.F_struc,mb.F_struc);
d = optimizer([F1,F3,x+y<=p],x+y,ops,p,[x;y]);
e = [b d];
me = e.model;
testCase.verifyEqual(me.K.s,[2 2 2]);

function test_polynomial_gradient_creation_order(testCase)
ops = testCase.TestData.options;
for ordering = 1:2
    x = sdpvar(1); f = x^3;
    y = sdpvar(1); f = f+y^3;
    if ordering == 2
        f = x*y^2+y^3+2*x^2;
        expectedValue = 4*9+27+2*16;
        expectedGradient = [9+16;24+27];
    else
        expectedValue = 64+27;
        expectedGradient = [48;27];
    end
    [~,~,~,m] = export([1<=x<=5,1<=y<=5],f,ops);
    m = yalmip2nonlinearsolver(m);
    global latest_xevaled latest_x_xevaled
    latest_xevaled = []; latest_x_xevaled = [];
    [v,g] = fmincon_fun([4;3],m);
    testCase.verifyEqual(v,expectedValue,'AbsTol',1e-12);
    testCase.verifyEqual(g(:),expectedGradient,'AbsTol',1e-12);
end

function test_sdp_presolve_retains_coefficients(testCase)
x = sdpvar(1); z = binvar(1);
m = prepareSDP([1 1;1 x]>=0,-x+z,x,testCase.TestData.options);
p = presolve_implied_from_SDP(m);
ix = find(m.used_variables==getvariables(x));
testCase.verifyEqual(p.lb(ix),10);
testCase.verifyEqual(p.F_struc,m.F_struc);

function test_sdp_presolve_quadratic_objective(testCase)
x = sdpvar(1); z = binvar(1);
m = prepareSDP([x 1;1 1]>=0,x^2+z,x,testCase.TestData.options);
p = presolve_implied_from_SDP(m);
ix = find(m.used_variables==getvariables(x));
testCase.verifyEqual([p.lb(ix) p.ub(ix)],[0 10]);

function m = prepareSDP(F,f,x,ops)
[~,~,~,m] = export([F,0<=x<=10],f,ops);
m = compile_nonlinear_table(m); m.KCut.f = []; m.KCut.l = [];
m = presolve_bounds_from_modelbounds(m,1); m.feasible = 1;
m.noninteger_variables = setdiff(1:length(m.c),[m.integer_variables m.binary_variables m.semicont_variables]);
