function tests = test_robust_14
tests = functiontests(localfunctions);

function test1(dummy)
sdpvar t x 
w = sdpvar(2,1);
M = momentmodel([w'*w == 1, w>=0.5])
P = sosmodel(sos(1 + w(1)*x + w(2)*x^2 - t),[],[],[t;w]);
sol = optimize([M,P,uncertain(M)],-t)
assert(sol.problem == 0);
assert(abs(value(t)-.625) <= 1e-3);

