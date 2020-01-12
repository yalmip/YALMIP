function tests = test_operator_optimizer2
tests = functiontests(localfunctions);

function test1(dummy)
% Tests issue #39

yalmip('clear')
sdpvar x
sdpvar z
sdpvar y
P = optimizer([-1 <= y <= 1,x <= 10+z*x],x^2,[],[y;z],x)
[~,err] = P{[0;5]};

assert(err == 0);

yalmip('clear')
sdpvar x
sdpvar y
sdpvar z
P = optimizer([-1 <= y <= 1,x <= 10+z*x],x^2,[],[y;z],x)
[~,err] = P{[0;5]};

assert(err == 0);
