function tests = test_bmibnb_nonlinearktkt
tests = functiontests(localfunctions);

function test1(testCase)

sdpvar x1 x2 y1 y2
x=[x1,x2];y=[y1,y2];
ops=sdpsettings('solver','bmibnb');

F = -x1^2-3*x2-4*y1+y2^2;
G = [x1^2+2*x2-4<=0;-x1<=0;-x2<=0];
f = 2*x1^2+y1^2-5*y2;
g = [-x1^2+2*x1-x2^2+2*y1-y2-3<=0;-x2-3*y1+4*y2+4<=0;-y1<=0;-y2<=0];
K = kkt(g,f,[x1 x2]);
sol = optimize([G,K,g],F);
testCase.assertTrue(sol.problem == 0)
testCase.assertTrue(abs(value(F)--12.6787) <= 1e-2) 