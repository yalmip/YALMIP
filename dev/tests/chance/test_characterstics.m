function tests = test_characterstics
tests = functiontests(localfunctions);

function test_exponential_single(testCase)
w = sdpvar(1);
sdpvar t x
gamma = sdpvar(1)

Model = [probability(3*w <= t) >= 1-gamma, 0 <= gamma <= 0.05
         uncertain(w,'exponential',3)];
     
optimize(Model,t,sdpsettings('debug',1,'fmincon.alg','sqp'))

% Confirm numerically
N = 1e6;
w1 = random('exponential',3,[1,N]);
wsample = w1;
estimated_probability = (nnz(3*wsample <= value(t)))/N

testCase.assertTrue(abs(value(t)-26.96)<=0.01 && abs(estimated_probability - 0.95) <= 0.01)


function test_exponential_multi(testCase)
w = sdpvar(2,1);
sdpvar t x
gamma = sdpvar(1)

Model = [probability([3 2]*w <= t) >= 1-gamma, 0 <= gamma <= 0.05
         uncertain(w,'exponential',[3;0.3])];
     
optimize(Model,t,sdpsettings('debug',1,'fmincon.alg','sqp'))

% Confirm numerically
N = 1e6;
w1 = random('exponential',3,[1,N]);
w2 = random('exponential',.3,[1,N]);
wsample = [w1;w2];
estimated_probability = (nnz([3 2]*wsample <= value(t)))/N

testCase.assertTrue(abs(value(t)-27.58)<=0.01 && abs(estimated_probability - 0.95) <= 0.01)



function test_logistic_single(testCase)
w = sdpvar(1);
sdpvar t x
gamma = sdpvar(1)

Model = [probability(2*w <= t) >= 1-gamma, 0 <= gamma <= 0.05
         uncertain(w,'logistic',0.1,2)];
     
optimize(Model,t,sdpsettings('debug',1,'fmincon.alg','sqp'))

% Confirm numerically
N = 1e6;
w1 = random('logistic',0.1,2,[1,N]);
wsample = w1;
estimated_probability = (nnz(2*wsample <= value(t)))/N;

testCase.assertTrue(abs(value(t)-11.9777)<=0.01 && abs(estimated_probability - 0.95) <= 0.01)


function test_logistic_multi(testCase)

w = sdpvar(2,1);
sdpvar t x

Model = [probability([2+x 3-x]*w <= t) >= 0.95,
         uncertain(w,'logistic',[0.1;-.3],[.05;3])];
     
optimize(Model,t,sdpsettings('debug',1,'fmincon.alg','sqp'))

% Confirm numerically
N = 1e6;
w1 = random('logistic',0.1,.05,[1,N]);
w2 = random('logistic',-.3,3,[1,N]);
wsample = [w1;w2];
xopt = value(x)
estimated_probability = (nnz([2+xopt 3-xopt]*wsample <= value(t)))/N;

testCase.assertTrue(abs(xopt-3)<=0.01 && abs(estimated_probability - 0.95) <= 0.01)


function test2(testCase)

w = sdpvar(2,1);
sdpvar t x 
gamma = sdpvar(1)

Model = [probability([2+x 3-x]*w <= t) >= 1-gamma,0 <=gamma <= 0.05,
         uncertain(w,'logistic',[0.1;-.3],[.05;3])];
     
optimize(Model,t,sdpsettings('debug',1,'fmincon.alg','sqp'))

% Confirm numerically
N = 1e6;
w1 = random('logistic',0.1,.05,[1,N]);
w2 = random('logistic',-.3,3,[1,N]);
wsample = [w1;w2];
xopt = value(x)
estimated_probability = (nnz([2+xopt 3-xopt]*wsample <= value(t)))/N;
abs(xopt-3)<=0.01
abs(estimated_probability - 0.95) <= 0.01
testCase.assertTrue(abs(xopt-3)<=0.01 && abs(estimated_probability - 0.95) <= 0.01)



function test3(testCase)

w = sdpvar(2,1);
sdpvar t x 

Model = [probability([2+x 3-x]*w <= t) >= 0.95,
         uncertain(w,'uniform',[1;-3],[2;3])];
     
optimize(Model,t,sdpsettings('debug',1,'fmincon.alg','sqp'))

% Confirm numerically
N = 1e4;
w1 = random('uniform',1,2,[1,N]);
w2 = random('uniform',-3,3,[1,N]);
wsample = [w1;w2];
xopt = value(x)
estimated_probability = (nnz([2+xopt 3-xopt]*wsample <= value(t)))/N

abs(xopt-3)<=0.01
abs(estimated_probability - 0.95) <= 0.01
testCase.assertTrue(abs(xopt-2.4769)<=0.01 && abs(estimated_probability - 0.95) <= 0.01)

