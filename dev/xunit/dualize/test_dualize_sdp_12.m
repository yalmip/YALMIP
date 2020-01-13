function tests = test_sdpvar_dualize_sdp_12
tests = functiontests(localfunctions);

function test1(dummy)

yalmip('clear')
N = 2;
X = sdpvar(N,N,'hermitian','complex');
x = sdpvar(3,1);
t = sdpvar;
obj = t;
F = [X>=0, X(1)+x(3) == 7];
F = [F,cone([1;t]),cone([30;X(:)])];
opts = sdpsettings('dualize',1);
sol1 = optimize(F,-obj,opts)
o1 = value(obj);
opts = sdpsettings('dualize',0);
sol2 = optimize(F,-obj,opts);
o2 = value(obj);
assert(sol1.problem == 0);
assert(sol2.problem == 0);
assert(abs(o1 - o2) <= 1e-4);
