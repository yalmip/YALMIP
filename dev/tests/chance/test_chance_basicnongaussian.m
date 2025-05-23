function tests = test_chance_basicnongaussian
tests = functiontests(localfunctions);

function test_scalar_unit_exponential(testCase)
yalmip('clear')
w = sdpvar(1,1);
sdpvar t
truth = 2.3078;

Model = [probability(2+3*w >= t) >= 0.95,uncertain(w,'exponential',2)];
optimize(Model,-t)
testCase.assertTrue(abs(value(t)-truth) <= 1e-3)
% nnz(2+3*random('exponential',2,1,1e6)>=value(t))/1e6

truth = -15.9744;
Model = [probability(2-3*w >= t) >= 0.95,uncertain(w,'exponential',2)];
optimize(Model,-t)
testCase.assertTrue(abs(value(t)-truth) <= 1e-3)
% nnz(2-3*random('exponential',2,1,1e6)>=value(t))/1e6

function test_scalar_unit_logistic(testCase)
yalmip('clear')
w = sdpvar(1,1);
sdpvar t
truth = -18.5;
Model = [probability(2+3*w >= t) >= 0.95,uncertain(w,'logistic',2,3)];
optimize(Model,-t)
testCase.assertTrue(abs(value(t)-truth) <= 1e-3)
%nnz(2+3*random('logistic',2,3,1,1e6)>=value(t))/1e6

truth = -30.5;
Model = [probability(2-3*w >= t) >= 0.95,uncertain(w,'logistic',2,3)];
optimize(Model,-t)
testCase.assertTrue(abs(value(t)-truth) <= 1e-3)
%nnz(2-3*random('logistic',2,3,1,1e6)>=value(t))/1e6

function test_different_scalar_uniform(testCase)
yalmip('clear')
w = sdpvar(1,1);
sdpvar t

truth = -9.1;
Model = [probability(2+3*w >= t) >= 0.95,uncertain(w,'uniform',-4,2)];
optimize(Model,-t)
testCase.assertTrue(abs(value(t)-truth) <= 1e-3)
%nnz(2+3*random('uniform',-4,2,1,1e6)>=value(t))/1e6

truth = -3.1;
Model = [probability(2-3*w >= t) >= 0.95,uncertain(w,'uniform',-4,2)];
optimize(Model,-t)
testCase.assertTrue(abs(value(t)-truth) <= 1e-3)
%nnz(2-3*random('uniform',-4,2,1,1e6)>=value(t))/1e6


function test_different_scalar_laplace(testCase)
yalmip('clear')
w = sdpvar(1,1);
mu = 2;
sigma = 3;
scale = sigma/sqrt(2);
sdpvar t

truth = -6.6536;
Model = [probability(2+3*w >= t) >= 0.95,uncertain(w,'laplace',mu,sigma)];
optimize(Model,-t)
testCase.assertTrue(abs(value(t)-truth) <= 1e-3)
%nnz(2+3*(mu+scale*(random('exponential',1,1,1e6)-random('exponential',1,1,1e6))) >= value(t))/1e6

truth = -18.6536;
Model = [probability(2-3*w >= t) >= 0.95,uncertain(w,'laplace',mu,sigma)];
optimize(Model,-t)
testCase.assertTrue(abs(value(t)-truth) <= 1e-3)
%nnz(2-3*(mu+scale*(random('exponential',1,1,1e6)-random('exponential',1,1,1e6))) >= value(t))/1e6


