function tests = test_misc_kktqp
tests = functiontests(localfunctions);

function test1(dummy)
Q = magic(5);
x = sdpvar(5,1);
OO = x'*Q*x;
optimize([-1 <= x <= 1],x'*Q*x,sdpsettings('solver','kktqp'))
assert(abs(value(OO)-(-80.9412))<=1e-2);
