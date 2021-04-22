function tests = test_operator_max
tests = functiontests(localfunctions);

function test1(dummy)
% Test for bug #345
a=sdpvar(1,1);
b=sdpvar(1,1);
C=[];
C=[C,a==-1];
C=[C,b==max(a,0)];
sol = optimize(C)
assert(abs(value(b)) <= 1e-4);

