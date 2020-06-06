function tests = test_bmibnb_specialfunctions
tests = functiontests(localfunctions);


%% erf functions
function test_erf(dummy)
sdpvar x
sol = optimize([-1 <= x <= 1], (erf(x-1/pi) - 1/pi)^2,sdpsettings('solver','bmibnb','bmibnb.uppersolver','fmincon'));
assert(sol.problem == 0)

function test_erfc(dummy)
sdpvar x
sol = optimize([-1 <= x <= 1], (erfc(x-1/pi) - 1/pi)^2,sdpsettings('solver','bmibnb','bmibnb.uppersolver','fmincon'));
assert(sol.problem == 0)

function test_erfcinv(dummy)
sdpvar x
sol = optimize([-1 <= x <= 1], (erfcinv(x-1/pi) - 1/pi)^2,sdpsettings('solver','bmibnb','bmibnb.uppersolver','fmincon'));
assert(sol.problem == 0)

function test_erfcx(dummy)
sdpvar x
sol = optimize([-1 <= x <= 1], (erfcx(x-1/pi) - 1/pi)^2,sdpsettings('solver','bmibnb','bmibnb.uppersolver','fmincon'));
assert(sol.problem == 0)

%% Gamma functions
function test_gamma(dummy)
sdpvar x
sol = optimize([-1 <= x <= 1], (gamma(x-1/pi) - 1/pi)^2,sdpsettings('solver','bmibnb','bmibnb.uppersolver','fmincon'));
assert(sol.problem == 0)

function test_gammaln(dummy)
sdpvar x
sol = optimize([-1 <= x <= 1], (gammaln(x-1/pi) - 1/pi)^2,sdpsettings('solver','bmibnb','bmibnb.uppersolver','fmincon'));
assert(sol.problem == 0)

function test_gammainc_x(dummy)
sdpvar x
sol = optimize([-1 <= x <= 1], (gammainc(x-1/pi,2) - 1/pi)^2,sdpsettings('solver','bmibnb','bmibnb.uppersolver','fmincon'));
assert(sol.problem == 0)

function test_gammainc_a(dummy)
sdpvar x
sol = optimize([-1 <= x <= 1], (gammainc(2,x-1/pi) - 1/pi)^2,sdpsettings('solver','bmibnb','bmibnb.uppersolver','fmincon'));
assert(sol.problem == 0)

function test_gammaincinv_x(dummy)
sdpvar x
sol = optimize([-1 <= x <= 1], (gammaincinv(0.5,x-1/pi) - 1/pi)^2,sdpsettings('solver','bmibnb','bmibnb.uppersolver','fmincon'));
assert(sol.problem == 0)

function testgammainv_a(dummy)
sdpvar x
sol = optimize([-1 <= x <= 1], (gammaincinv(x-1/pi,0.5) - 1/pi)^2,sdpsettings('solver','bmibnb','bmibnb.uppersolver','fmincon'));
assert(sol.problem == 0)

%% Cos functions
function test_cos(dummy)
sdpvar x
sol = optimize([-1 <= x <= 1], (cos(x-1/pi) - 1/pi)^2,sdpsettings('solver','bmibnb','bmibnb.uppersolver','fmincon'));
assert(sol.problem == 0)

function test_acos(dummy)
sdpvar x
sol = optimize([-1 <= x <= 1], (acos(x-1/pi) - 1/pi)^2,sdpsettings('solver','bmibnb','bmibnb.uppersolver','fmincon'));
assert(sol.problem == 0)

function test_cosh(dummy)
sdpvar x
sol = optimize([-1 <= x <= 1], (cosh(x-1/pi) - 1/pi)^2,sdpsettings('solver','bmibnb','bmibnb.uppersolver','fmincon'));
assert(sol.problem == 0)

function test_acosh(dummy)
sdpvar x
sol = optimize([-1 <= x <= 1], (acosh(x-1/pi) - 1/pi)^2,sdpsettings('solver','bmibnb','bmibnb.uppersolver','fmincon'));
assert(sol.problem == 0)

%% sin functions
function test_sin(dummy)
sdpvar x
sol = optimize([-1 <= x <= 1], (sin(x-1/pi) - 1/pi)^2,sdpsettings('solver','bmibnb','bmibnb.uppersolver','fmincon'));
assert(sol.problem == 0)

function test_asin(dummy)
sdpvar x
sol = optimize([-1 <= x <= 1], (asin(x-1/pi) - 1/pi)^2,sdpsettings('solver','bmibnb','bmibnb.uppersolver','fmincon'));
assert(sol.problem == 0)

function test_sinh(dummy)
sdpvar x
sol = optimize([-1 <= x <= 1], (sinh(x-1/pi) - 1/pi)^2,sdpsettings('solver','bmibnb','bmibnb.uppersolver','fmincon'));
assert(sol.problem == 0)

function test_asinh(dummy)
sdpvar x
sol = optimize([-1 <= x <= 1], (asinh(x-1/pi) - 1/pi)^2,sdpsettings('solver','bmibnb','bmibnb.uppersolver','fmincon'));
assert(sol.problem == 0)

%% tan functions
function test_tan(dummy)
sdpvar x
sol = optimize([-1 <= x <= 1], (tan(x-1/pi) - 1/pi)^2,sdpsettings('solver','bmibnb','bmibnb.uppersolver','fmincon'));
assert(sol.problem == 0)

function test_atan(dummy)
sdpvar x
sol = optimize([-1 <= x <= 1], (atan(x-1/pi) - 1/pi)^2,sdpsettings('solver','bmibnb','bmibnb.uppersolver','fmincon'));
assert(sol.problem == 0)

function test_tanh(dummy)
sdpvar x
sol = optimize([-1 <= x <= 1], (tanh(x-1/pi) - 1/pi)^2,sdpsettings('solver','bmibnb','bmibnb.uppersolver','fmincon'));
assert(sol.problem == 0)

function test_atanh(dummy)
sdpvar x
sol = optimize([-1 <= x <= 1], (atanh(x-1/pi) - 1/pi)^2,sdpsettings('solver','bmibnb','bmibnb.uppersolver','fmincon'));
assert(sol.problem == 0)

%% sec functions
function test_sec(dummy)
sdpvar x
sol = optimize([-1 <= x <= 1], (sec(x-1/pi) - 1/pi)^2,sdpsettings('solver','bmibnb','bmibnb.uppersolver','fmincon'));
assert(sol.problem == 0)

function test_asec(dummy)
sdpvar x
sol = optimize([-1 <= x <= 1], (asec(x-1/pi) - 1/pi)^2,sdpsettings('solver','bmibnb','bmibnb.uppersolver','fmincon'));
assert(sol.problem == 0)

function test_sech(dummy)
sdpvar x
sol = optimize([-1 <= x <= 1], (sech(x-1/pi) - 1/pi)^2,sdpsettings('solver','bmibnb','bmibnb.uppersolver','fmincon'));
assert(sol.problem == 0)

function test_asech(dummy)
sdpvar x
sol = optimize([-1 <= x <= 1], (asech(x-1/pi) - 1/pi)^2,sdpsettings('solver','bmibnb','bmibnb.uppersolver','fmincon'));
assert(sol.problem == 0)

%% cot functions
function test_cot(dummy)
sdpvar x
sol = optimize([-1 <= x <= 1], (cot(x-1/pi) - 1/pi)^2,sdpsettings('solver','bmibnb','bmibnb.uppersolver','fmincon'));
assert(sol.problem == 0)

function test_acot(dummy)
sdpvar x
sol = optimize([-1 <= x <= 1], (acot(x-1/pi) - 1/pi)^2,sdpsettings('solver','bmibnb','bmibnb.uppersolver','fmincon'));
assert(sol.problem == 0)

function test_coth(dummy)
sdpvar x
sol = optimize([-1 <= x <= 1], (coth(x-1/pi) - 1/pi)^2,sdpsettings('solver','bmibnb','bmibnb.uppersolver','fmincon'));
assert(sol.problem == 0)

function test_acoth(dummy)
sdpvar x
sol = optimize([-1 <= x <= 2], (acoth(x-1/pi) - 1/pi)^2,sdpsettings('solver','bmibnb','bmibnb.uppersolver','fmincon'));
assert(sol.problem == 0)

%% csc functions
function test_csc(dummy)
% Use this to test sign-switching singularity effect on bounds etc
sdpvar x
sol = optimize([-1 <= x <= 1], (csc(x-1/pi) - 1/pi)^2,sdpsettings('solver','bmibnb','bmibnb.uppersolver','fmincon'));
assert(sol.problem == 0)

function test_acsc(dummy)
% Disconnected domain!
%sdpvar x
%sol = optimize([-1 <= x <= 1], (acsc(x-1/pi) - 1/pi)^2,sdpsettings('solver','bmibnb','bmibnb.uppersolver','fmincon'));
%assert(sol.problem == 0)

function test_csch(dummy)
sdpvar x
sol = optimize([-1 <= x <= 1], (csch(x-1/pi) - 1/pi)^2,sdpsettings('solver','bmibnb','bmibnb.uppersolver','fmincon'));
assert(sol.problem == 0)

function test_acsch(dummy)
sdpvar x
sol = optimize([-1 <= x <= 1], (acsch(x-1/pi) - 1/pi)^2,sdpsettings('solver','bmibnb','bmibnb.uppersolver','fmincon'));
assert(sol.problem == 0)

function test_airy(dummy)
sdpvar x
sol = optimize([-1 <= x <= 1], (airy(x-1/pi) - 1/pi)^2,sdpsettings('solver','bmibnb','bmibnb.uppersolver','fmincon'));
assert(sol.problem == 0)

function test_beta(dummy)
sdpvar x
sol = optimize([-1 <= x <= 1], (beta(x-1/pi,.5) - 1/pi)^2,sdpsettings('solver','bmibnb','bmibnb.uppersolver','fmincon'));
assert(sol.problem == 0)

function test_expint(dummy)
sdpvar x
sol = optimize([-1 <= x <= 1], (expint(x-1/pi) - 1/pi)^2,sdpsettings('solver','bmibnb','bmibnb.uppersolver','fmincon'));
assert(sol.problem == 0)

function test_sqrtm(dummy)
sdpvar x
sol = optimize([-1 <= x <= 1], (sqrtm(x-1/pi) - 1/pi)^2,sdpsettings('solver','bmibnb','bmibnb.uppersolver','fmincon'));
assert(sol.problem == 0)

