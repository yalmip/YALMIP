function tests = test_rcone
tests = functiontests(localfunctions);

function test1(dummy)
x = sdpvar(1);
y = sdpvar(1);
z = sdpvar(3,1);

optimize([z>=1, rcone(z,x,y)],x+y);
assert(abs(value(x*y)-1.5) <= 1e-4);