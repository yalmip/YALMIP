function tests = test_global_ex14_1_2
tests = functiontests(localfunctions);

function test1(dummy)

yalmip('clear')
sdpvar x1 x2 x3 x4 x5 x6 objvar

F = ([]);
F = F + ( - x6 + objvar == 0); 
F = F + ( x1*x2 + x1 - 3*x5 == 0);
F = F + ( 2.8845e-6*sqr(x2) + 4.4975e-7*x2 + 2*x1*x2 + x1 + 0.000545176668613029*x2* x3 + 3.40735417883143e-5*x2*x4 + x2*sqr(x3) - 10*x5 - x6 <= 0);
F = F + ( (-2.8845e-6*sqr(x2)) - 4.4975e-7*x2 - 2*x1*x2 - x1 - 0.000545176668613029* x2*x3 - 3.40735417883143e-5*x2*x4 - x2*sqr(x3) + 10*x5 - x6 <= 0);
F = F + ( 0.386*sqr(x3) + 0.000410621754172864*x3 + 0.000545176668613029*x2*x3 + 2* x2*sqr(x3) - 8*x5 - x6 <= 0); 
F = F + ( (-0.386*sqr(x3)) - 0.000410621754172864*x3 - 0.000545176668613029*x2*x3 - 2*x2*sqr(x3) + 8*x5 - x6 <= 0); 
F = F + ( 2*sqr(x4) + 3.40735417883143e-5*x2*x4 - 40*x5 - x6 <= 0);
F = F + ( (-2*sqr(x4)) - 3.40735417883143e-5*x2*x4 + 40*x5 - x6 <= 0);
F = F + ( 9.615e-7*sqr(x2) + 4.4975e-7*x2 + 0.193*sqr(x3) + 0.000410621754172864*x3 + sqr(x4) + x1*x2 + x1 + 0.000545176668613029*x2*x3 + 3.40735417883143e-5 *x2*x4 + x2*sqr(x3) - x6 <= 1);
F = F + ( (-9.615e-7*sqr(x2)) - 4.4975e-7*x2 - 0.193*sqr(x3) - 0.000410621754172864 *x3 - sqr(x4) - x1*x2 - x1 - 0.000545176668613029*x2*x3 - 3.40735417883143e-5*x2*x4 - x2*sqr(x3) - x6 <= -1);
F = F + (0.0001 <= [x1 x2 x3 x4 x5] <= 100);
F = F + (0.0001*ones(1,5) <= [x1 x2 x3 x4 x5] <= 100);

sol = optimize(F,objvar,sdpsettings('solver','bmibnb'));
assert(sol.problem == 0)
assert(value(objvar) <= 1e-2)

function y = sqr(x)
y = x*x;