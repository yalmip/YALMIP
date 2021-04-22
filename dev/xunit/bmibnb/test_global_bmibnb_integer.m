function tests = test_global_bmibnb_integer
tests = functiontests(localfunctions);

function test1(dummy)

clear sin % To avoid aby problems from code above
sdpvar x
strange = @(x) sdpfun(x,'@(x) sin(10*x)+abs(sin(x))+x');
obj = strange(x);
sol=optimize((integer(x)) + (-pi <= x <= pi),strange(x),sdpsettings('solver','bmibnb'));

assert(sol.problem==0)
assert(abs(value(x)--2) <= 1e-4)
