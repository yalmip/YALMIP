function tests = test_chance_sampler
tests = functiontests(localfunctions);

function test1(testCase)
yalmip('clear')
w = sdpvar(2,1);
Model = [uncertain(w,'exponential',[2;1])];
W = sample(w);
testCase.assertTrue(isequal(size(W),[2 1]))

function test2(testCase)
yalmip('clear')
w = sdpvar(2,1);
Model = [uncertain(w,'normal',[1;1],[2;2])];
W = sample(w,3);
testCase.assertTrue(isequal(size(W),[2 3]))

function test3(testCase)
yalmip('clear')
w = sdpvar(2,1);
Model = [uncertain(w,'normal mixture',{[1;1],[-10;-5]},{[1;1],[2;2]},[.25 .75;0.99 0.01])];
W = sample(w,10);
testCase.assertTrue(isequal(size(W),[2 10]))

function test4(testCase)
yalmip('clear')
w = sdpvar(2,1);
Model = [uncertain(w(1),'normal mixture',{[1],[-10]},{[1],[2]},[.25 .75]),
         uncertain(w(2),'normal mixture',{[1],[-5]},{[1],[2]},[.9 .1])];
W = sample(w,10);
testCase.assertTrue(isequal(size(W),[2 10]))

function test5(testCase)
yalmip('clear')
w = sdpvar(2,1);
u = sdpvar(3,1);
Model = [uncertain(w(1),'normal mixture',{[1],[-10]},{[1],[2]},[.25 .75]),
         uncertain(w(2),'normal mixture',{[1],[-5]},{[1],[2]},[.9 .1]);
         uncertain(u,'exponential mixture',{[1],[5]},[.9 .1])];
W = sample(w,10);
U = sample(u,10);
testCase.assertTrue(isequal(size(W),[2 10]))
testCase.assertTrue(isequal(size(U),[3 10]))

function test6(testCase)
yalmip('clear')
w = sdpvar(2,1);
u = sdpvar(3,1);
Model = [uncertain(w(1),'normal mixture',{[1],[-10]},{[1],[2]},[.25 .75]),
         uncertain(w(2),'normal mixture',{[1],[-5]},{[1],[2]},[.9 .1]);
         uncertain(u,'exponential mixture',{[1],[5]},[.9 .1])];
A = sample([w;u(1)-u(2)],10);
testCase.assertTrue(isequal(size(A),[3 10]))

function test7(testCase)
yalmip('clear');
w1 = sdpvar(3,1);
w2 = sdpvar(3,1);
w3 = sdpvar(3,1);

Model = [uncertain(w1,'mvnrnd mixture',{[1;2;3],[5;6;7]},{1e-4*eye(3), 1e-4*eye(3) + 1e-4*ones(3)},[.75 .25])]
Model = [uncertain(w3,'normal mixture',{[1000;1000;1000],[-1000;-1000;-1000]},{1e-4*[1;1;1],1e-4*[1;1;1]},[.5 .5])]
Model = [uncertain(w2,'mvnrnd mixture',{[8;9;10],[11;12;13]},{1e-4*eye(3), 1e-4*eye(3) + 1e-4*ones(3)},[.05 .95])]
A = sample([w1;w3;w2]);
testCase.assertTrue(isequal(size(A),[9 1]))



