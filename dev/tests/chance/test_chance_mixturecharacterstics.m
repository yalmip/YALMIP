function tests = test_chance_mixturecharacterstics
tests = functiontests(localfunctions);

function test_logistic_multi_mixture(testCase)

w = sdpvar(2,1);
sdpvar t x
gamma = sdpvar(1)

mu = {[0.1;-0.3], [0.5;0.7]};
theta = {[0.05;3], [0.04;2]};
weights = [0.25, 0.75];

Model = [probability([2+x 3-x]*w <= t) >= 1-gamma, 
         0 <= gamma <= 0.05, 
         uncertain(w,'logistic mixture',mu,theta,weights)];
     
% This should call the characterstic framework     
optimize(Model,t,sdpsettings('debug',1,'fmincon.alg','sqp'))

% Confirm numerically
N = 1e6;
% First element
w11 = random('logistic',mu{1}(1),theta{1}(1),[1,weights(1)*N]);
w12 = random('logistic',mu{2}(1),theta{2}(1),[1,weights(2)*N]);
w1 = [w11 w12];
% Second element
w21 = random('logistic',mu{1}(2),theta{1}(2),[1,weights(1)*N]);
w22 = random('logistic',mu{2}(2),theta{2}(2),[1,weights(2)*N]);
w2 = [w21 w22];
% Two components independantly drawn
wsample = [w1;w2];
wsample = [w1(randperm(length(wsample)));w2(randperm(length(wsample)))];

xopt = value(x)
estimated_probability = (nnz([2+xopt 3-xopt]*wsample <= value(t)))/length(wsample)

testCase.assertTrue(abs(xopt-3)<=0.01 && abs(estimated_probability - 0.95) <= 0.01)


function test_gaussian_mixture(testCase)

yalmip('clear');
gamma = sdpvar(1);
sdpvar u w x
a = 1/4;
b = 1/4;
c = 1/4;
d = 1/4;
P1 = probability((1+x)*w<=u);
truth = 7.1581;
Model = [uncertain(w,'normal mixture',{-3,0,3,2},{10,2,1,2},{a, b, c, d}),P1 >= 1-gamma, gamma <= 0.04, x >= 0]
optimize(Model,u,sdpsettings('chance.char','yes','debug',1,'fmincon.alg','sqp'))

testCase.assertTrue(abs(value(u)-truth)<=0.01)



