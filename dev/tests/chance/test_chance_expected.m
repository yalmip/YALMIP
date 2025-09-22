function tests = test_chance_expected
tests = functiontests(localfunctions);

function test1(testCase)
yalmip('clear')
w = sdpvar(1,1);
Model = [uncertain(w,'normal',2,3)];
E = expected(w);
testCase.assertTrue(E == 2);

function test2(testCase)
yalmip('clear')
w = sdpvar(2,1);
Model = [uncertain(w,'normal',[2;3],[4;5])];
E = expected(w);
testCase.assertTrue(all(E == [2;3]));

function test3(testCase)
yalmip('clear')
w1 = sdpvar(2,1);
w2 = sdpvar(2,1);
Model = [uncertain(w1,'normal',[4;5],[4;5])];
         uncertain(w2,'normal',[2;3],[4;5])
E = expected(w2);
testCase.assertTrue(all(E == [2;3]));

function test4(testCase)
yalmip('clear')
w1 = sdpvar(2,1);
w2 = sdpvar(2,1);
Model = [uncertain(w1,'normal',[4;5],[4;5])];
         uncertain(w2,'normal',[2;3],[4;5])
E = expected([3 4]*w2);
testCase.assertTrue(all(E == [3 4]*[2;3]));

function test5(testCase)
yalmip('clear')
sdpvar x
w1 = sdpvar(2,1);
w2 = sdpvar(2,1);
Model = [uncertain(w1,'normal',[4;5],[4;5])];
         uncertain(w2,'normal',[2;3],[4;5])
E = expected(x + [3+x 4*x]*w2);
assign(x,2);
testCase.assertTrue(all(value(E) == 2 + [3+2 4*2]*[2;3]));