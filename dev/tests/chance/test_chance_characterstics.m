function tests = test_characterstics
tests = functiontests(localfunctions);

function test_exponential_single(testCase)
w = sdpvar(1);
sdpvar t x
gamma = sdpvar(1)

Model = [probability(3*w <= t) >= 1-gamma, 0 <= gamma <= 0.05, uncertain(w,'exponential',3)];
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

Model = [probability([3 2]*w <= t) >= 1-gamma, 0 <= gamma <= 0.05, uncertain(w,'exponential',[3;0.3])];
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

Model = [probability(2*w <= t) >= 1-gamma, 0 <= gamma <= 0.05, uncertain(w,'logistic',0.1,2)];
% Force use of characteristic functions
optimize(Model,t,sdpsettings('debug',1,'solver','fmincon','fmincon.alg','sqp','chance.characteristic','yes'))

% Confirm numerically
N = 1e6;
w1 = random('logistic',0.1,2,[1,N]);
wsample = w1;
estimated_probability = (nnz(2*wsample <= value(t)))/N;

testCase.assertTrue(abs(value(t)-11.9777)<=0.01 && abs(estimated_probability - 0.95) <= 0.01)

function test_cauchy_single(testCase)
yalmip('clear')
w=sdpvar(1,1);
sdpvar t g
truth = -16.9413;
Model = [probability(w >= t) >= 1-g,0 <= g <= 0.05,uncertain(w,'cauchy',2,3)];
optimize(Model,-t,sdpsettings('debug',1,'solver','fmincon','fmincon.alg','sqp','chance.characteristic','yes'))
testCase.assertTrue(abs(value(t)-truth) <= 1e-3)
%d = makedist('tLocationScale','mu',2,'sigma',3,'nu',1)
%r = random(d,10000,1);
%nnz(r >= value(t))/length(r)

function test_laplace_single(testCase)
yalmip('clear')
w=sdpvar(1,1);
sdpvar t g
truth = -3.2562;
Model = [probability(w >= t) >= 1-g,0 <= g <= 0.05,uncertain(w,'laplace',0,2)];
optimize(Model,-t,sdpsettings('debug',1,'solver','fmincon','fmincon.alg','','chance.characteristic','yes'))
testCase.assertTrue(abs(value(t)-truth) <= 1e-3)
if 0
    r = random('exponential',2/sqrt(2),1,100000);
    r = r - random('exponential',2/sqrt(2),1,100000);
    nnz(r >= value(t))/length(r)
end

function test_laplace_multi(testCase)
yalmip('clear')
w=sdpvar(2,1);
sdpvar t g
truth = -13.8154;
Model = [probability(2*w(1)+3*w(2) >= t) >= 1-g,0 <= g <= 0.05,uncertain(w,'laplace',[1;2],[3;4])];
optimize(Model,-t,sdpsettings('debug',1,'solver','fmincon','fmincon.alg','','chance.characteristic','yes'))
testCase.assertTrue(abs(value(t)-truth) <= 1e-3)
if 0
    r1 = random('exponential',3/sqrt(2),1,100000);
    r1 = r1 - random('exponential',3/sqrt(2),1,100000) + 1;
    r2 = random('exponential',4/sqrt(2),1,100000);
    r2 = r2 - random('exponential',4/sqrt(2),1,100000) + 2;
    nnz(2*r1+3*r2 >= value(t))/length(r)
end


function test_cauchy_multi(testCase)
yalmip('clear')
w=sdpvar(2,1);
sdpvar t g
truth = -84.7063;
Model = [probability([2 3]*w >= t) >= 1-g,0 <= g <= 0.05,uncertain(w,'cauchy',2,3)];
optimize(Model,-t,sdpsettings('debug',1,'solver','fmincon','fmincon.alg','sqp','chance.characteristic','yes'))
testCase.assertTrue(abs(value(t)-truth) <= 1e-3)
%d = makedist('tLocationScale','mu',2,'sigma',3,'nu',1)
%r = random(d,2,1e5);
%nnz([2 3]*r >= value(t))/length(r)

function test_cauchy_multi_2(testCase)
yalmip('clear')
w=sdpvar(2,1);
sdpvar t g 
x = sdpvar(100,1);
truth = -84.7063;
Model = [probability([2+sum(x) 3-sum(x)]*w >= t) >= 1-g,0 <= g <= 0.05,uncertain(w,'cauchy',2,3)];
optimize(Model,-t,sdpsettings('debug',1,'solver','fmincon','fmincon.alg','sqp','chance.characteristic','yes'))
testCase.assertTrue(abs(value(t)-truth) <= 1e-3)
%d = makedist('tLocationScale','mu',2,'sigma',3,'nu',1)
%r = random(d,2,1e5);
%nnz([2 3]*r >= value(t))/length(r)


function test_logistic_multi(testCase)

w = sdpvar(2,1);
sdpvar t x
gamma = sdpvar(1)

Model = [probability([2+x 3-x]*w <= t) >= 1-gamma, 0 <= gamma <= 0.05, uncertain(w,'logistic',[0.1;-.3],[.05;3])];
     
% This should call the characterstic framework     
optimize(Model,t,sdpsettings('debug',1,'fmincon.alg','sqp'))

% Confirm numerically
N = 1e6;
w1 = random('logistic',0.1,.05,[1,N]);
w2 = random('logistic',-.3,3,[1,N]);
wsample = [w1;w2];
xopt = value(x)
estimated_probability = (nnz([2+xopt 3-xopt]*wsample <= value(t)))/N;
testCase.assertTrue(abs(xopt-3)<=0.01 && abs(estimated_probability - 0.95) <= 0.01)

function test_uniform_multi(testCase)

w = sdpvar(2,1);
sdpvar t x 
gamma = sdpvar(1)

% This example struggles with SQP solver if gamma is a decision variable
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
testCase.assertTrue(abs(xopt-2.4769)<=0.01 && abs(estimated_probability - 0.95) <= 0.01)

