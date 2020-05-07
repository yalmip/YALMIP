function tests = test_misc_binmodel
tests = functiontests(localfunctions);

function test1(dummy)
sdpvar y
binvar x

[xy,model] = binmodel(x*y,[2 <= y <= 5]);
optimize(model,-xy)
assert(abs(value(xy)-5) <= 1e-4);
[xy,model] = binmodel(x*y,abs(y) <= 5);
optimize(model,xy)
assert(abs(value(xy)-(-5))<=1e-4);
