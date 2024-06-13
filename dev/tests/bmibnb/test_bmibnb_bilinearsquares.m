function tests = test_bmibnb_bilinearsquares
tests = functiontests(localfunctions);

function test1(testCase)

sdpvar x1 x2
Model = [-3 <= [x1 x2] <= 3];
obj = x1-7*x2+x1^2*x2^2-.2*x1^4*x2^4+.1*x1^6*x2^6;
sol = optimize(Model,obj,sdpsettings('solver','bmibnb','savesolveroutput',1,'debug',1,'bmibnb.branchinor',0,'verbose',0));
testCase.assertTrue(abs(value(obj)--21.0279) <= 1e-3)
testCase.assertTrue(sol.solveroutput.nodes <= 20)

function test2(testCase)
sdpvar x1 x2 x3
Model = [-3 <= [x1 x2] <= 3,x3==x1*x2 ];
obj = x1-7*x2+x3^2-.2*x3^4+.1*x3^6;
sol = optimize(Model,obj,sdpsettings('solver','bmibnb','savesolveroutput',1,'debug',1,'verbose',0));

testCase.assertTrue(abs(value(obj)--21.0279) <= 1e-3)
testCase.assertTrue(sol.solveroutput.nodes <= 20)