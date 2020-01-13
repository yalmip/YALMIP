function tests = test_geomean
tests = functiontests(localfunctions);

function test1(dummy)
% Test real vector geomean, length == 2^n
randn('seed',1234);
rand('seed',1234);
A = randn(15,2);
b = rand(15,1)*5;
x = sdpvar(2,1);
obj = geomean(b-A*x);
optimize([],-obj);
assert(norm(value(x') - [-0.05519469470525   0.26970610928222]) <= 1e-3);
assert(norm(value(obj) - 1.83896843735621) <= 1e-3);

% Test real vector geomean, length == 2^n
randn('seed',1234);
rand('seed',1234);
A = randn(16,2);
b = rand(16,1)*5;
x = sdpvar(2,1);
obj = geomean(b-A*x);
optimize([],-obj);
assert(norm(value(x') - [ -0.01148934254297  -0.20720944929269]) <= 1e-3);
assert(abs(value(obj) - 1.93924577959868) <= 1e-3);

% Test real vector geomean, length == 1
randn('seed',1234);
rand('seed',1234);
A = randn(1,2);
b = rand(1,1)*5;
x = sdpvar(2,1);
obj = geomean(b-A*x);
sol = optimize([],-obj);
assert(sol.problem==2);

% Test real matrix geomean, length ~2^n
randn('seed',1234);
rand('seed',1234);
D = randn(5,5);
P = sdpvar(5,5);
obj = geomean(P);
optimize((P <= D*D'),-obj);
assert(abs(value(obj) - 2.00333629658259) <= 1e-4);
 
% Test real matrix geomean, length == 2^n
randn('seed',1234);
rand('seed',1234);
D = randn(8,8);
P = sdpvar(8,8);
obj = geomean(P);
optimize((P <= D*D'),-obj);
assert(abs(value(obj) - 3.32199302165511) <= 1e-4);

% Test real matrix geomean, length == 2
randn('seed',1234);
rand('seed',1234);
D = randn(2,2);
P = sdpvar(2,2);
obj = geomean(P);
optimize((P <= D*D'),-obj);
assert(abs(value(obj) - 2.02896175488410) <= 1e-4);