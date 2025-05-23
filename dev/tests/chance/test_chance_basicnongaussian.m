function tests = test_chance_basicnongaussian
tests = functiontests(localfunctions);

function test_scalar_unit_exponential(testCase)
yalmip('clear')
w=sdpvar(1,1);
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
w=sdpvar(1,1);
sdpvar t
truth = 2.3078;


Model = [probability(w >= t) >= 0.95,uncertain(w,'logistic',1,1)];
optimize(Model,-t)
testCase.assertTrue(abs(value(t)-truth) <= 1e-3)
% nnz(random('logistic',1,1,1,1e6)>=value(t))/1e6

Model = [probability(2+3*w >= t) >= 0.95,uncertain(w,'logistic',2,3)];
optimize(Model,-t)
testCase.assertTrue(abs(value(t)-truth) <= 1e-3)
% nnz(2+3*random('logistic',2,3,1,1e6)>=value(t))/1e6

truth = -15.9744;
Model = [probability(2-3*w >= t) >= 0.95,uncertain(w,'exponential',2,3)];
optimize(Model,-t)
testCase.assertTrue(abs(value(t)-truth) <= 1e-3)
% nnz(2-3*random('exponential',2,3,1,1e6)>=value(t))/1e6



function test_different_scalar_gaussian(testCase)
yalmip('clear')
w=sdpvar(1,1);
sdpvar t
 
% Test that std.dev vs variance is correctly used
Model = [probability(w >= t) >= 0.95,uncertain(w,'normal',0,2)];
optimize(Model,-t)
testCase.assertTrue(abs(value(t)--3.2897) <= 1e-3)

Model = [probability(w >= t) >= 0.95,uncertain(w,'mvnrnd',0,2^2)];
optimize(Model,-t)
testCase.assertTrue(abs(value(t)--3.2897) <= 1e-3)

Model = [probability(w >= t) >= 0.95,uncertain(w,'mvnrndfactor',0,2)];
optimize(Model,-t)
testCase.assertTrue(abs(value(t)--3.2897) <= 1e-3)

function test_scalarmultivariate_gaussian(testCase)
yalmip('clear')
w=sdpvar(2,1);
c = [2;3];
sdpvar t
mu = [1;2];
sigma = [3;4];
truth = -14.068;

% nnz(c'*(mu+sigma.*random('normal',0,1,2,1e6))>=value(t))/1e6
Model = [probability(c'*w >= t) >= 0.95,uncertain(w,'normal',mu,sigma)];
optimize(Model,-t)
testCase.assertTrue(abs(value(t)-truth) <= 1e-3)

Model = [probability(c'*w >= t) >= 0.95,uncertain(w,'mvnrnd',mu,diag(sigma.^2))];
optimize(Model,-t)
testCase.assertTrue(abs(value(t)-truth) <= 1e-3)

Model = [probability(c'*w >= t) >= 0.95,uncertain(w,'mvnrndfactor',mu,diag(sigma))];
optimize(Model,-t)
testCase.assertTrue(abs(value(t)-truth) <= 1e-3)

Model = [probability(c'*w >= t) >= 0.95,uncertain(w,'mvnrnd',mu,sigma.^2)];
optimize(Model,-t)
testCase.assertTrue(abs(value(t)-truth) <= 1e-3)

Model = [probability(c'*w >= t) >= 0.95,uncertain(w,'mvnrndfactor',mu,sigma)];
optimize(Model,-t)
testCase.assertTrue(abs(value(t)-truth) <= 1e-3)

function test_trulymultivariate_gaussian(testCase)
yalmip('clear')
w=sdpvar(2,1);
c = [2;3];
sdpvar t
mu = [1;2];

% nnz(c'*(mu+R*random('normal',0,1,2,1e6))>=value(t))/1e6
truth = -9.485;
R = [2 1;1 2];C = R'*R;
Model = [probability(c'*w >= t) >= 0.95,uncertain(w,'mvnrndfactor',mu,R)];
optimize(Model,-t)
testCase.assertTrue(abs(value(t)-truth) <= 1e-3)

Model = [probability(c'*w >= t) >= 0.95,uncertain(w,'mvnrnd',mu,C)];
optimize(Model,-t)
testCase.assertTrue(abs(value(t)-truth) <= 1e-3)

Model = [probability(c'*(mu + R*w) >= t) >= 0.95,uncertain(w,'mvnrnd',[0;0],eye(2))];
optimize(Model,-t)
testCase.assertTrue(abs(value(t)-truth) <= 1e-3)

Rv = sdpvar(2,2,'full');
Model = [probability(c'*w >= t) >= 0.95,uncertain(w,'mvnrndfactor',mu,Rv), R == Rv];
optimize(Model,-t)
testCase.assertTrue(abs(value(t)-truth) <= 1e-3)

function test_concatenate_scalar_gaussian(testCase)
yalmip('clear')
sdpvar w1 w2 t
sdpvar t

Model = [probability(w1 >= t) >= 0.5,
         probability(w2 >= t) >= 0.5, 
         uncertain(w1,'normal',0,eye(1)),
         uncertain(w2,'normal',0,eye(1))];
optimize(Model,-t)
testCase.assertTrue(abs(value(t)) <= 1e-3)
