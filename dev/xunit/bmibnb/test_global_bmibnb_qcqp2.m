function tests = test_global_bmibnb_qcqp2
tests = functiontests(localfunctions);

function test1(dummy)

x1 = sdpvar(1,1);
x2 = sdpvar(1,1);
x3 = sdpvar(1,1);
x4 = sdpvar(1,1);
x5 = sdpvar(1,1);
x6 = sdpvar(1,1);

obj = -25*(x1-2)^2-(x2-2)^2-(x3-1)^2-(x4-4)^2-(x5-1)^2-(x6-4)^2;

F = ((x3-3)^2+x4>=4) + ((x5-3)^2+x6>=4);
F = F + (x1-3*x2<=2) + (-x1+x2<=2) + (x1+x2>=2);
F = F + (6>=x1+x2>=2);
F = F + (1<=x3<=5) + (0<=x4<=6)+(1<=x5<=5)+(0<=x6<=10)+(x1>=0)+(x2>=0);
sol = optimize(F,obj,sdpsettings('solver','bmibnb'))

assert(abs(value(obj)--310) <= 1e-6)