function tests = test_sdpvar_hessian
tests = functiontests(localfunctions);

function test1(testCase)

sdpvar x y
testCase.assertTrue(all(all(zeros(2)==hessian(x,[x y]))));
testCase.assertTrue(all(all(zeros(2)==hessian(x+y,[x y]))));
testCase.assertTrue(all(all([0 1;1 0]==hessian(x*y,[x y]))));
testCase.assertTrue(all(all([2 0;0 0]==hessian(x^2,[x y]))));
testCase.assertTrue(all(all([2]==hessian(x^2+y^2,[x]))));
testCase.assertTrue(all(all([2]==hessian(x^2+y^2,[y]))));
testCase.assertTrue(all(all([0]==hessian(y^2,[x]))));






 
