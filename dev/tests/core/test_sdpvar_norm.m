function tests = test_sdpvar_norm
tests = functiontests(localfunctions);

function test1(testCase)
% To improve performance, we don't introduce auxilliary variables for
% elements which are fixed
sdpvar x y
assign(x,1);
testCase.assertTrue(abs(value(norm([x;-3],1)) - 4)<= 1e-8)
optimize(norm([x;-3],1)<=y,y);
testCase.assertTrue(norm(value(y) - 3)<= 1e-8)
testCase.assertTrue(norm(value(x) - 0)<= 1e-8)
