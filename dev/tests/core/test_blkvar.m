function tests = test_blkvar
tests = functiontests(localfunctions);

function test1(testCase)

x = sdpvar(1);
Yblk = blkvar;
Yblk(1,1) = x;
Yblk(1,2) = 2;
Yblk(2,1) = 2;
Yblk(2,2) = 1+x;
Y = [x 2;2 1+x];
testCase.assertTrue(isequal(Y - Yblk,zeros(2)));
testCase.assertTrue(isequal(Y+Y'-(Yblk+Yblk'),zeros(2)));

Yblk = blkvar;
Yblk(1,1) = x;
Yblk(2,2) = 1+x;
Yblk(3,3) = 1-x;
testCase.assertTrue(isequal(diag([x 1+x 1-x]) - Yblk,zeros(3)));


