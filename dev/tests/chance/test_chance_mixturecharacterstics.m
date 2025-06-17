function tests = test_chance_mixturecharacterstics
tests = functiontests(localfunctions);

function test_logistic_multi(testCase)

w = sdpvar(2,1);
sdpvar t x
gamma = sdpvar(1)

mu = {[0.1;-0.3], [0.5;0.7]};
theta = {[0.05;3], [0.04;2]};
weights = {0.5, 0.5};

Model = [probability([2+x 3-x]*w <= t) >= 1-gamma, 
         0 <= gamma <= 0.05, 
         uncertain(w,'logistic mixture',mu,theta,weights)];
     
% This should call the characterstic framework     
optimize(Model,t,sdpsettings('debug',1,'fmincon.alg','sqp'))

% Confirm numerically
N = 1e6;
% equal probability of drawing from first or second mix component for both
% elements, so totally for combinations of equal likelihood
w11 = random('logistic',mu{1}(1),theta{1}(1),[1,N]);
w12 = random('logistic',mu{2}(1),theta{2}(1),[1,N]);
w21 = random('logistic',mu{1}(2),theta{1}(2),[1,N]);
w22 = random('logistic',mu{2}(2),theta{2}(2),[1,N]);
wsample = [ [w11;w12] [w21;w12] [w11;w22] [w12;w22]];

xopt = value(x)
estimated_probability = (nnz([2+xopt 3-xopt]*wsample <= value(t)))/(4*N);
testCase.assertTrue(abs(xopt-3)<=0.01 && abs(estimated_probability - 0.95) <= 0.01)
