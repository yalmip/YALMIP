function tests = test_global_bmibnb_gamsrobot5
tests = functiontests(localfunctions);

function test1(dummy)

x1 = sdpvar(1,1);
x2 = sdpvar(1,1);
x3 = sdpvar(1,1);
x4 = sdpvar(1,1);
x5 = sdpvar(1,1);
x6 = sdpvar(1,1);
x7 = sdpvar(1,1);
x8 = sdpvar(1,1);
x = [x1;x2;x3;x4;x5;x6;x7;x8];
F = (-1<=x<=1);
F = F + ( 0.004731*x1*x3 - 0.1238*x1 - 0.3578*x2*x3 - 0.001637*x2 - 0.9338*x4 + x7 == 0.3571);
F = F + (0.2238*x1*x3 + 0.2638*x1 + 0.7623*x2*x3 - 0.07745*x2 - 0.6734*x4 - x7 == 0.6022);
F = F + ( x6*x8 + 0.3578*x1 + 0.004731*x2 == 0);
F = F + ( - 0.7623*x1 + 0.2238*x2 == -0.3461);
F = F + (x1^2 + x2^2 == 1);
F = F + (x3^2 + x4^2 == 1);
F = F + (x5^2 + x6^2 == 1);
F = F + (x7^2 + x8^2 == 1);
obj = sum(x);

sol = optimize(F,obj,sdpsettings('solver','bmibnb'))

assert(abs(value(obj)--3.5295) <= 1e-4)
