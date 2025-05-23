function tests = test_chance_basicnongaussian
tests = functiontests(localfunctions);

function test_scalar_unit_exponential(testCase)
yalmip('clear')
w = sdpvar(1,1);
sdpvar t s
lambda = 2;
truth = 2.3078;

Model = [probability(2+3*w >= t) >= 0.95,uncertain(w,'exponential',lambda)];
optimize(Model,-t)
testCase.assertTrue(abs(value(t)-truth) <= 1e-3)
% nnz(2+3*random('exponential',lambda,1,1e6)>=value(t))/1e6

Model = [probability(2+3*w >= t) >= 1-s, s <= 0.05, uncertain(w,'exponential',lambda)];
optimize(Model,-t)
testCase.assertTrue(abs(value(t)-truth) <= 1e-3)

truth = -15.9744;
Model = [probability(2-3*w >= t) >= 0.95,uncertain(w,'exponential',lambda)];
optimize(Model,-t)
testCase.assertTrue(abs(value(t)-truth) <= 1e-3)
% nnz(2-3*random('exponential',lambda,1,1e6)>=value(t))/1e6

Model = [probability(2-3*w >= t) >= 1-s, s <= 0.05, uncertain(w,'exponential',lambda)];
optimize(Model,-t)
testCase.assertTrue(abs(value(t)-truth) <= 1e-3)

function test_scalar_unit_weibull(testCase)
yalmip('clear')
w = sdpvar(1,1);
sdpvar t s
a = 1.8;
b = 2;
truth = 3.223;

Model = [probability(2+3*w >= t) >= 0.95,uncertain(w,'weibull',a,b)];
optimize(Model,-t)
testCase.assertTrue(abs(value(t)-truth) <= 1e-3)
% nnz(2+3*random('weibull',a,b,1,1e6)>=value(t))/1e6

Model = [probability(2+3*w >= t) >= 1-s, s <= 0.05, uncertain(w,'weibull',a,b)];
optimize(Model,-t)
testCase.assertTrue(abs(value(t)-truth) <= 1e-3)

truth = -7.3464;
Model = [probability(2-3*w >= t) >= 0.95,uncertain(w,'weibull',a,b)];
optimize(Model,-t)
testCase.assertTrue(abs(value(t)-truth) <= 1e-3)
% nnz(2-3*random('weibull',a,b,1,1e6)>=value(t))/1e6

Model = [probability(2-3*w >= t) >= 1-s, s <= 0.05, uncertain(w,'weibull',a,b)];
optimize(Model,-t)
testCase.assertTrue(abs(value(t)-truth) <= 1e-3)

function test_scalar_unit_gamma(testCase)
yalmip('clear')
w = sdpvar(1,1);
sdpvar t s
a = 1.8;
b = 2;
truth = 3.6696;

Model = [probability(2+3*w >= t) >= 0.95,uncertain(w,'gamma',a,b)];
optimize(Model,-t)
testCase.assertTrue(abs(value(t)-truth) <= 1e-3)
% nnz(2+3*random('gamma',a,b,1,1e6)>=value(t))/1e6

Model = [probability(2+3*w >= t) >= 1-s, s <= 0.05, uncertain(w,'gamma',a,b)];
optimize(Model,-t)
testCase.assertTrue(abs(value(t)-truth) <= 1e-3)

truth = -24.495;
Model = [probability(2-3*w >= t) >= 0.95,uncertain(w,'gamma',a,b)];
optimize(Model,-t)
testCase.assertTrue(abs(value(t)-truth) <= 1e-3)
% nnz(2-3*random('gamma',a,b,1,1e6)>=value(t))/1e6

Model = [probability(2-3*w >= t) >= 1-s, s <= 0.05, uncertain(w,'gamma',a,b)];
optimize(Model,-t)
testCase.assertTrue(abs(value(t)-truth) <= 1e-3)


function test_scalar_unit_logistic(testCase)
yalmip('clear')
w = sdpvar(1,1);
sdpvar t s
truth = -18.5;
Model = [probability(2+3*w >= t) >= 0.95,uncertain(w,'logistic',2,3)];
optimize(Model,-t)
testCase.assertTrue(abs(value(t)-truth) <= 1e-3)
%nnz(2+3*random('logistic',2,3,1,1e6)>=value(t))/1e6

truth = -18.5;
Model = [probability(2+3*w >= t) >= 1-s, s <= 0.05, uncertain(w,'logistic',2,3)];
optimize(Model,-t)
testCase.assertTrue(abs(value(t)-truth) <= 1e-3)

truth = -30.5;
Model = [probability(2-3*w >= t) >= 0.95,uncertain(w,'logistic',2,3)];
optimize(Model,-t)
testCase.assertTrue(abs(value(t)-truth) <= 1e-3)
%nnz(2-3*random('logistic',2,3,1,1e6)>=value(t))/1e6

truth = -30.5;
Model = [probability(2-3*w >= t) >= 1-s, s <= 0.05, uncertain(w,'logistic',2,3)];
optimize(Model,-t)
testCase.assertTrue(abs(value(t)-truth) <= 1e-3)


function test_different_scalar_uniform(testCase)
yalmip('clear')
w = sdpvar(1,1);
sdpvar t s

truth = -9.1;
Model = [probability(2+3*w >= t) >= 0.95,uncertain(w,'uniform',-4,2)];
optimize(Model,-t)
testCase.assertTrue(abs(value(t)-truth) <= 1e-3)
%nnz(2+3*random('uniform',-4,2,1,1e6)>=value(t))/1e6

Model = [probability(2+3*w >= t) >= 1-s, s <= 0.05,uncertain(w,'uniform',-4,2)];
optimize(Model,-t)
testCase.assertTrue(abs(value(t)-truth) <= 1e-3)

truth = -3.1;
Model = [probability(2-3*w >= t) >= 0.95,uncertain(w,'uniform',-4,2)];
optimize(Model,-t)
testCase.assertTrue(abs(value(t)-truth) <= 1e-3)
%nnz(2-3*random('uniform',-4,2,1,1e6)>=value(t))/1e6

Model = [probability(2-3*w >= t) >= 1-s, s <= 0.05,uncertain(w,'uniform',-4,2)];
optimize(Model,-t)
testCase.assertTrue(abs(value(t)-truth) <= 1e-3)

function test_different_scalar_laplace(testCase)
yalmip('clear')
w = sdpvar(1,1);
mu = 2;
sigma = 3;
scale = sigma/sqrt(2);
sdpvar t s

truth = -6.6536;
Model = [probability(2+3*w >= t) >= 0.95,uncertain(w,'laplace',mu,sigma)];
optimize(Model,-t)
testCase.assertTrue(abs(value(t)-truth) <= 1e-3)
%nnz(2+3*(mu+scale*(random('exponential',1,1,1e6)-random('exponential',1,1,1e6))) >= value(t))/1e6

Model = [probability(2+3*w >= t) >= 1-s, s <= 0.05, uncertain(w,'laplace',mu,sigma)];
optimize(Model,-t)
testCase.assertTrue(abs(value(t)-truth) <= 1e-3)

truth = -18.6536;
Model = [probability(2-3*w >= t) >= 0.95,uncertain(w,'laplace',mu,sigma)];
optimize(Model,-t)
testCase.assertTrue(abs(value(t)-truth) <= 1e-3)
%nnz(2-3*(mu+scale*(random('exponential',1,1,1e6)-random('exponential',1,1,1e6))) >= value(t))/1e6

Model = [probability(2-3*w >= t) >= 1-s, s <= 0.05, uncertain(w,'laplace',mu,sigma)];
optimize(Model,-t)
testCase.assertTrue(abs(value(t)-truth) <= 1e-3)

function test_different_scalar_studentt(testCase)
yalmip('clear')
w = sdpvar(1,1);
nu = 2;
sdpvar t s

truth = -6.76;
Model = [probability(2+3*w >= t) >= 0.95,uncertain(w,'t',nu)];
optimize(Model,-t)
testCase.assertTrue(abs(value(t)-truth) <= 1e-3)
%nnz(2+3*random('t',nu,1,1e6) >= value(t))/1e6

Model = [probability(2+3*w >= t) >= 1-s, s <= 0.05, uncertain(w,'t',nu)];
optimize(Model,-t)
testCase.assertTrue(abs(value(t)-truth) <= 1e-3)

function test_different_scalar_tlocationscale(testCase)
yalmip('clear')
w = sdpvar(1,1);
sdpvar t s

truth = -18.1477;
Model = [probability(2+3*w >= t) >= 0.95,uncertain(w,'tlocationScale',3,5,6)];
optimize(Model,-t)
testCase.assertTrue(abs(value(t)-truth) <= 1e-3)
%nnz(2+3*random('tlocationScale',3,5,6,1,1e6) >= value(t))/1e6

Model = [probability(2+3*w >= t) >= 1-s, s <= 0.05, uncertain(w,'tlocationScale',3,5,6)];
optimize(Model,-t)
testCase.assertTrue(abs(value(t)-truth) <= 1e-3)

function test_different_scalar_stable_symmetric(testCase)
yalmip('clear')
w = sdpvar(1,1);
a = 2;
b = 0; % b = 0 leads to symmetric
c = 1;
d = 0;
sdpvar t s

truth = -4.9785;
Model = [probability(2+3*w >= t) >= 0.95,uncertain(w,'stable',a,b,c,d)];
optimize(Model,-t)
testCase.assertTrue(abs(value(t)-truth) <= 1e-3)
%nnz(2+3*random('stable',a,b,c,d,1,1e6) >= value(t))/1e6

Model = [probability(2+3*w >= t) >= 1-s, s <= 0.05, uncertain(w,'stable',a,b,c,d)];
optimize(Model,-t)
testCase.assertTrue(abs(value(t)-truth) <= 1e-3)
