function tests = test_bmibnb_specialfunctions
tests = functiontests(localfunctions);

%% erf functions
function test_erf(testCase)
sdpvar x
sol = optimize([-1 <= x <= 1], (erf(x-1/pi) - 1/pi)^2,sdpsettings('solver','bmibnb','bmibnb.uppersolver','fmincon'));
testCase.assertTrue(sol.problem == 0)

function test_erfc(testCase)
sdpvar x
sol = optimize([-1 <= x <= 1], (erfc(x-1/pi) - 1/pi)^2,sdpsettings('solver','bmibnb','bmibnb.uppersolver','fmincon'));
testCase.assertTrue(sol.problem == 0)

function test_erfcinv(testCase)
sdpvar x
sol = optimize([-1 <= x <= 1], (erfcinv(x-1/pi) - 1/pi)^2,sdpsettings('solver','bmibnb','bmibnb.uppersolver','fmincon'));
testCase.assertTrue(sol.problem == 0)

function test_erfcx(testCase)
sdpvar x
sol = optimize([-1 <= x <= 1], (erfcx(x-1/pi) - 1/pi)^2,sdpsettings('solver','bmibnb','bmibnb.uppersolver','fmincon'));
testCase.assertTrue(sol.problem == 0)

function test_expint(testCase)
sdpvar x
sol = optimize([-1 <= x <= 1], (expint(x-1/pi) - 1/pi)^2,sdpsettings('solver','bmibnb','bmibnb.uppersolver','fmincon'));
testCase.assertTrue(sol.problem == 0)

%% Gamma functions
function test_gamma(testCase)
sdpvar x
sol = optimize([-1 <= x <= 1], (gamma(x-1/pi) - 1/pi)^2,sdpsettings('solver','bmibnb','bmibnb.uppersolver','fmincon'));
testCase.assertTrue(sol.problem == 0)

function test_gammaln(testCase)
sdpvar x
sol = optimize([-1 <= x <= 1], (gammaln(x-1/pi) - 1/pi)^2,sdpsettings('solver','bmibnb','bmibnb.uppersolver','fmincon'));
testCase.assertTrue(sol.problem == 0)

function test_gammainc_x(testCase)
sdpvar x
sol = optimize([-1 <= x <= 1], (gammainc(x-1/pi,2) - 1/pi)^2,sdpsettings('solver','bmibnb','bmibnb.uppersolver','fmincon'));
testCase.assertTrue(sol.problem == 0)

function test_gammainc_a(testCase)
sdpvar x
sol = optimize([-1 <= x <= 1], (gammainc(2,x-1/pi) - 1/pi)^2,sdpsettings('solver','bmibnb','bmibnb.uppersolver','fmincon'));
testCase.assertTrue(sol.problem == 0)

function test_gammaincinv_x(testCase)
sdpvar x
sol = optimize([-1 <= x <= 1], (gammaincinv(0.5,x-1/pi) - 1/pi)^2,sdpsettings('solver','bmibnb','bmibnb.uppersolver','fmincon'));
testCase.assertTrue(sol.problem == 0)

function testgammainv_a(testCase)
sdpvar x
sol = optimize([-1 <= x <= 1], (gammaincinv(x-1/pi,0.5) - 1/pi)^2,sdpsettings('solver','bmibnb','bmibnb.uppersolver','fmincon'));
testCase.assertTrue(sol.problem == 0)

%% Cos functions
function test_cos(testCase)
sdpvar x
sol = optimize([-1 <= x <= 1], (cos(x-1/pi) - 1/pi)^2,sdpsettings('solver','bmibnb','bmibnb.uppersolver','fmincon'));
testCase.assertTrue(sol.problem == 0)

function test_acos(testCase)
sdpvar x
sol = optimize([-1 <= x <= 1], (acos(x-1/pi) - 1/pi)^2,sdpsettings('solver','bmibnb','bmibnb.uppersolver','fmincon'));
testCase.assertTrue(sol.problem == 0)

function test_cosh(testCase)
sdpvar x
sol = optimize([-1 <= x <= 1], (cosh(x-1/pi) - 1/pi)^2,sdpsettings('solver','bmibnb','bmibnb.uppersolver','fmincon'));
testCase.assertTrue(sol.problem == 0)

function test_acosh(testCase)
sdpvar x
sol = optimize([-5 <= x <= 5], (acosh(x-1/pi) - 1/pi)^2,sdpsettings('solver','bmibnb','bmibnb.uppersolver','fmincon'));
testCase.assertTrue(sol.problem == 0)

%% sin functions
function test_sin(testCase)
sdpvar x
sol = optimize([-1 <= x <= 1], (sin(x-1/pi) - 1/pi)^2,sdpsettings('solver','bmibnb','bmibnb.uppersolver','fmincon'));
testCase.assertTrue(sol.problem == 0)

function test_asin(testCase)
sdpvar x
sol = optimize([-1 <= x <= 1], (asin(x-1/pi) - 1/pi)^2,sdpsettings('solver','bmibnb','bmibnb.uppersolver','fmincon'));
testCase.assertTrue(sol.problem == 0)

function test_sinh(testCase)
sdpvar x
sol = optimize([-1 <= x <= 1], (sinh(x-1/pi) - 1/pi)^2,sdpsettings('solver','bmibnb','bmibnb.uppersolver','fmincon'));
testCase.assertTrue(sol.problem == 0)

function test_asinh(testCase)
sdpvar x
sol = optimize([-1 <= x <= 1], (asinh(x-1/pi) - 1/pi)^2,sdpsettings('solver','bmibnb','bmibnb.uppersolver','fmincon'));
testCase.assertTrue(sol.problem == 0)

%% tan functions
function test_tan(testCase)
sdpvar x
sol = optimize([-1 <= x <= 1], (tan(x-1/pi) - 1/pi)^2,sdpsettings('solver','bmibnb','bmibnb.uppersolver','fmincon'));
testCase.assertTrue(sol.problem == 0)

function test_atan(testCase)
sdpvar x
sol = optimize([-1 <= x <= 1], (atan(x-1/pi) - 1/pi)^2,sdpsettings('solver','bmibnb','bmibnb.uppersolver','fmincon'));
testCase.assertTrue(sol.problem == 0)

function test_tanh(testCase)
sdpvar x
sol = optimize([-1 <= x <= 1], (tanh(x-1/pi) - 1/pi)^2,sdpsettings('solver','bmibnb','bmibnb.uppersolver','fmincon'));
testCase.assertTrue(sol.problem == 0)

function test_atanh(testCase)
sdpvar x
sol = optimize([-1 <= x <= 1], (atanh(x-1/pi) - 1/pi)^2,sdpsettings('solver','bmibnb','bmibnb.uppersolver','fmincon'));
testCase.assertTrue(sol.problem == 0)

%% sec functions
function test_sec(testCase)
sdpvar x
sol = optimize([-1 <= x <= 1], (sec(x-1/pi) - 1/pi)^2,sdpsettings('solver','bmibnb','bmibnb.uppersolver','fmincon'));
testCase.assertTrue(sol.problem == 0)

function test_asec(testCase)
sdpvar x
sol = optimize([-1 <= x <= 1], (asec(x-1/pi) - 1/pi)^2,sdpsettings('solver','bmibnb','bmibnb.uppersolver','fmincon'));
testCase.assertTrue(sol.problem == 0)

function test_sech(testCase)
sdpvar x
sol = optimize([-1 <= x <= 1], (sech(x-1/pi) - 1/pi)^2,sdpsettings('solver','bmibnb','bmibnb.uppersolver','fmincon'));
testCase.assertTrue(sol.problem == 0)

function test_asech(testCase)
sdpvar x
sol = optimize([-1 <= x <= 1], (asech(x-1/pi) - 1/pi)^2,sdpsettings('solver','bmibnb','bmibnb.uppersolver','fmincon'));
testCase.assertTrue(sol.problem == 0)

%% cot functions
function test_cot(testCase)
sdpvar x
sol = optimize([-1 <= x <= 1], (cot(x-1/pi) - 1/pi)^2,sdpsettings('solver','bmibnb','bmibnb.uppersolver','fmincon'));
testCase.assertTrue(sol.problem == 0)

function test_acot(testCase)
sdpvar x
sol = optimize([-1 <= x <= 1], (acot(x-1/pi) - 1/pi)^2,sdpsettings('solver','bmibnb','bmibnb.uppersolver','fmincon'));
testCase.assertTrue(sol.problem == 0)

function test_coth(testCase)
sdpvar x
sol = optimize([-1 <= x <= 1], (coth(x-1/pi) - 1/pi)^2,sdpsettings('solver','bmibnb','bmibnb.uppersolver','fmincon'));
testCase.assertTrue(sol.problem == 0)

function test_acoth(testCase)
sdpvar x
sol = optimize([-1 <= x <= 2], (acoth(x-1/pi) - 1/pi)^2,sdpsettings('solver','bmibnb','bmibnb.uppersolver','fmincon'));
testCase.assertTrue(sol.problem == 0)

%% csc functions
function test_csc(testCase)
% Use this to test sign-switching singularity effect on bounds etc
sdpvar x
sol = optimize([-1 <= x <= 1], (csc(x-1/pi) - 1/pi)^2,sdpsettings('solver','bmibnb','bmibnb.uppersolver','fmincon'));
testCase.assertTrue(sol.problem == 0)

function test_acsc(testCase)
sdpvar x
sol = optimize([-10 <= x <= 10], (acsc(x-1/pi) - 1/pi)^2,sdpsettings('solver','bmibnb','bmibnb.uppersolver','fmincon'));
testCase.assertTrue(sol.problem == 0)

function test_csch(testCase)
sdpvar x
sol = optimize([-1 <= x <= 1], (csch(x-1/pi) - 1/pi)^2,sdpsettings('solver','bmibnb','bmibnb.uppersolver','fmincon'));
testCase.assertTrue(sol.problem == 0)

function test_acsch(testCase)
sdpvar x
sol = optimize([-1 <= x <= 1], (acsch(x-1/pi) - 1/pi)^2,sdpsettings('solver','bmibnb','bmibnb.uppersolver','fmincon'));
testCase.assertTrue(sol.problem == 0)

function test_airy(testCase)
sdpvar x
sol = optimize([-1 <= x <= 1], (airy(x-1/pi) - 1/pi)^2,sdpsettings('solver','bmibnb','bmibnb.uppersolver','fmincon'));
testCase.assertTrue(sol.problem == 0)

function test_beta(testCase)
sdpvar x
sol = optimize([-1 <= x <= 10], (beta(x-1/pi,.5) - 2)^2,sdpsettings('solver','bmibnb','bmibnb.uppersolver','fmincon'));
testCase.assertTrue(sol.problem == 0)

sol = optimize([-1 <= x <= 10], (beta(2,x-1/pi) - 2)^2,sdpsettings('solver','bmibnb','bmibnb.uppersolver','fmincon'));
testCase.assertTrue(sol.problem == 0)

function test_sqrtm(testCase)
sdpvar x
sol = optimize([-1 <= x <= 1], (sqrtm(x-1/pi) - 1/pi)^2,sdpsettings('solver','bmibnb','bmibnb.uppersolver','fmincon'));
testCase.assertTrue(sol.problem == 0)

