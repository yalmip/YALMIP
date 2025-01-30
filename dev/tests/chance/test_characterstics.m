function tests = test_characterstics
tests = functiontests(localfunctions);

function test1(testCase)

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
testCase.assertTrue(abs(xopt-3)<=0.01 && abs(estimated_probability - 0.95) <= 0.01)

