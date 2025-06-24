function tests = test_chance_basicjointgaussian
tests = functiontests(localfunctions);

function test_gaussian_joint_3x3(testCase)
yalmip('clear')
sdpvar x1 x2 x3 w1 w2 w3
truth = 1.83319
Model = [probability([w1 <= x1+x2+x3, 
                      w2 <= x1+x2+x3,
                      w3 <= x1+x2+x3],'joint') >= 1-0.1,
         uncertain(w1,'normal',0,1);   
         uncertain(w2,'normal',0,1)
         uncertain(w3,'normal',0,1)];
optimize(Model,x1+x2+x3)
testCase.assertTrue(value(x1+x2+x3)-truth <= 1e-3)

