function tests = test_bnb_lpqp
tests = functiontests(localfunctions);

function test_milp(testCase)
rand('seed',1234);

a = [1 2 3 4 5 6]';
t = (0:0.02:2*pi)';
x = [sin(t) sin(2*t) sin(3*t) sin(4*t) sin(5*t) sin(6*t)];
y = x*a+(-4+8*rand(length(x),1));

a_hat = intvar(6,1);

residuals = y-x*a_hat;
bound = sdpvar(length(residuals),1);
F = (-bound <= residuals <= bound);
ops = sdpsettings('solver','bnb','verbose',0);

% Test QP
obj = sum(bound);
sol = optimize(F,obj,ops);
testCase.assertTrue(sol.problem == 0);
testCase.assertTrue(abs(value(obj) - 6.168422746718130e+002) <= 1e-5);

function test_miqp(testCase)
rand('seed',1234);

a = [1 2 3 4 5 6]';
t = (0:0.02:2*pi)';
x = [sin(t) sin(2*t) sin(3*t) sin(4*t) sin(5*t) sin(6*t)];
y = x*a+(-4+8*rand(length(x),1));

a_hat = intvar(6,1);

residuals = y-x*a_hat;
bound = sdpvar(length(residuals),1);
F = (-bound <= residuals <= bound);
ops = sdpsettings('solver','bnb','verbose',0);

% Test QP
obj = residuals'*residuals;
sol = optimize((residuals <= 50),obj,ops);
testCase.assertTrue(sol.problem == 0);
testCase.assertTrue(abs(value(obj) - 1.605058709613011e+003) <= 1e-5);

function test_migp(testCase)
x = sdpvar(7,1);

% Data
a     = ones(7,1);
alpha = ones(7,1);
beta  = ones(7,1);
gamma = ones(7,1);
f = [1 0.8 1 0.7 0.7 0.5 0.5]';
e = [1 2 1 1.5 1.5 1 2]';
Cout6 = 10;
Cout7 = 10;

% Model
C = alpha+beta.*x;
A = sum(a.*x);
P = sum(f.*e.*x);
R = gamma./x;

D1 = R(1)*(C(4));
D2 = R(2)*(C(4)+C(5));
D3 = R(3)*(C(5)+C(7));
D4 = R(4)*(C(6)+C(7));
D5 = R(5)*(C(7));
D6 = R(6)*Cout6;
D7 = R(7)*Cout7;

% Constraints
F = (x >= 1) + (P <= 20) + (A <= 100);

% Objective
D = max([(D1+D4+D6),(D1+D4+D7),(D2+D4+D6),(D2+D4+D7),(D2+D5+D7),(D3+D5+D6),(D3+D7)]);

% Solve integer problem
ops = sdpsettings('solver','bnb','verbose',0);
sol = optimize(F+(integer(x)),D,ops);

testCase.assertTrue(sol.problem == 0);
testCase.assertTrue(all(abs(value(x) - [ 2     3     3     3     2     3     3]') <= 1e-3));
testCase.assertTrue(abs(value(D)-(8+1/3)) <= 1e-3);

function test_presolve(testCase)

intvar x y z
sol = optimize([x==0,cone([x;y-1;y+1])],z,sdpsettings('solver','bnb'));
testCase.assertTrue(sol.problem == 1);

sol = optimize([x==0,cone([z;y-1;y+1])],z,sdpsettings('solver','bnb'));
testCase.assertTrue(sol.problem == 0);
testCase.assertTrue(abs(value(z)-2) <= 1e-3);

sol = optimize([x==0,[x y;y x]>=0],x,sdpsettings('solver','bnb'));
testCase.assertTrue(sol.problem == 0);
testCase.assertTrue(abs(value(x)-0) <= 1e-3);

sol = optimize([x==0,[x y;y x]>=0,cone([y;z-1;z+1])],x,sdpsettings('solver','bnb'));
testCase.assertTrue(sol.problem == 1);
