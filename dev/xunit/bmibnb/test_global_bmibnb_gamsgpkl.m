function tests = test_global_bmibnb_gamsgpkl
tests = functiontests(localfunctions);

function test1(dummy)
x1 = sdpvar(1,1);
x2 = sdpvar(1,1);
t = sdpvar(1,1);

F = ([x1;x2]>=0);
F = F + (- x1 + x2 <= 1);
F = F + (x1 - x2 <= 1);
F = F + (- x1 + 2*x2 <= 3);
F = F + (2*x1 - x2 <= 3);
obj = (2*x1 - 2*x1*x1 + 2*x1*x2 + 3*x2 - 2*x2*x2);
F = F + (obj<=t);
sol = optimize(F,obj,sdpsettings('solver','bmibnb'))
assert(abs(value(obj)--3) <= 1e-4)