function tests = test_bnb_milp_1
tests = functiontests(localfunctions);

function test1(dummy)
rand('seed',1234);

a = [1 2 3 4 5 6]';
t = (0:0.02:2*pi)';
x = [sin(t) sin(2*t) sin(3*t) sin(4*t) sin(5*t) sin(6*t)];
y = x*a+(-4+8*rand(length(x),1));

a_hat = intvar(6,1);

residuals = y-x*a_hat;
bound = sdpvar(length(residuals),1);
F = (-bound <= residuals <= bound);
ops = sdpsettings('solver','bnb');

% Test QP
obj = sum(bound);
sol = optimize(F,obj,ops);
assert(sol.problem == 0);
assert(abs(value(obj) - 6.168422746718130e+002) <= 1e-5);

