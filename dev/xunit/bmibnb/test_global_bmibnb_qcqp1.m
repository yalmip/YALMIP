function tests = test_global_bmibnb_qcqp1
tests = functiontests(localfunctions);

function test1(dummy)

x1 = sdpvar(1,1);
x2 = sdpvar(1,1);
x3 = sdpvar(1,1);

obj = -2*x1+x2-x3;

F = (x1*(4*x1-4*x2+4*x3-20)+x2*(2*x2-2*x3+9)+x3*(2*x3-13)+24>=0);
F = F + (4-(x1+x2+x3)>=0);
F = F + (6-(3*x2+x3)>=0);
F = F + (x1>=0);
F = F + (2-x1>=0);
F = F + (x2>=0);
F = F + (x3>=0);
F = F + (3-x3>=0);

sol = optimize(F,obj,sdpsettings('solver','bmibnb'))

assert(abs(value(obj)--4) <= 1e-3)