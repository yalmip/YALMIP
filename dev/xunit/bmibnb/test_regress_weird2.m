function tests = test_regress_weird2
tests = functiontests(localfunctions);

function test1(dummy)

clear sin 
sdpvar x
strange = @(x) sdpfun(x,'@(x) sin(10*x)+abs(sin(x))+x');
sol = optimize((-pi <= x <= pi),strange(x),sdpsettings('solver','bmibnb'));

assert(sol.problem == 0)
assert(abs(value(strange(x))--3.21) <= 1e-2) 
clear sin
