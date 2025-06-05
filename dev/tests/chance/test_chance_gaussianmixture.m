function tests = test_chance_gaussianmixture
tests = functiontests(localfunctions);


function test_gaussian_scalar(testCase)
% Trivial case
yalmip('clear');
gamma = sdpvar(1);
sdpvar u w
a = 1/4;
b = 3/4;
P1 = probability(w<=u);
Model = [uncertain(w,'normal mixture',{-5,5},{1,sqrt(3)},{a, b}),P1 >= 1-gamma, gamma <= 0.04]
optimize(Model,u)
testCase.assertTrue(abs(value(u)-7.7945) <= 1e-3)

% Numerical verification
if 0
    M = 1e6;
    x = -5+randn(1,a*M);
    y = 5+sqrt(3)*randn(1,b*M);
    W = [x y];
    [value(gamma) 1-nnz(W < value(u))/M]
end


function test_gaussian_jointconstraint(testCase)

% Bonferroni on constraints, mixture done by YALMIP
yalmip('clear');
gamma = sdpvar(1);
sdpvar gamma1 gamma2 u w
a = 1/4;
b = 3/4;
P1 = probability(8-u-w>=0);
P2 = probability(8+u+w>=0);
gamma = gamma1 + gamma2;
Model = [uncertain(w,'normal mixture',{-5,5},{1,sqrt(3)},{a, b}),
    P1 >= 1-gamma1,
    P2 >= 1-gamma2]
optimize(Model,gamma)
testCase.assertTrue(abs(value(gamma)-0.0135) <= 1e-3)

% Numerical verification
if 0
    M = 1e6;
    x = -5+randn(1,a*M);
    y = 5+sqrt(3)*randn(1,b*M);
    W = [x y];
    % Check
    [value(gamma1) nnz(value(u)+W > 8)/M]
    [value(gamma2) nnz(value(u)+W < -8)/M]
    [value(gamma1+gamma2) nnz(value(u)+W < -8 | value(u)+W > 8)/M]
end



function test_gaussian_multiterms(testCase)

% Case with two elementwise gaussian mixtures each with two components,
% meaning the linear component will be expanded to 2*2 = 4 components
yalmip('clear');
gamma = sdpvar(1);
u = sdpvar(1);
w = sdpvar(2,1);
a = 0.25;
b = 0.75;
P1 = probability(w(1)+w(2)<=u);
Model = [uncertain(w,'normal mixture',{[-5;-3],[5;2]},{[1;2],[sqrt(3);3]},{a, b}),
    P1 >= 1-gamma, gamma <= 0.04]
optimize(Model,u)
testCase.assertTrue(abs(value(u)-12.0844) <= 1e-3)

% Numerical verification
if 0
    N = 1e6;
    W11 = -5 + 1.*randn(1,N);
    W12 = 5 + sqrt(3).*randn(1,N);
    W21 = -3 + 2.*randn(1,N);
    W22 = 2 + 3.*randn(1,N);
    W = [W11;W21];
    i = find(rand(1,N)>a);
    W(1,i) = W12(i);
    i = find(rand(1,N)>a);
    W(2,i) = W22(i);
    Y = [1 1]*W;
    % Check
    [value(gamma) 1-nnz(Y < value(u))/length(Y)]
    
    % Mixture sum manual from scratch
    sdpvar  w11 w12 w21 w22 gamma1 gamma2 gamma3 gamma4
    P1 = probability(w11+w21<=u);
    P2 = probability(w12+w21<=u);
    P3 = probability(w11+w22<=u);
    P4 = probability(w12+w22<=u);
    Model = [uncertain(w11,'normal',-5, 1),
        uncertain(w12,'normal',5, sqrt(3)),
        uncertain(w21,'normal',-3,2),
        uncertain(w22,'normal',2, 3),
        %a^2*P1+a*b*P2+b*a*P3+b^2*P4 >= 1-gamma, gamma <= 0.04]
        P1 >= 1-gamma1,
        P2 >= 1-gamma2,
        P3 >= 1-gamma3,
        P4 >= 1-gamma4,
        %a^2*P1+a*b*P2+a*b*P3+b^2*P4 >= 1-gamma,
        gamma <= 0.04, gamma == a^2*gamma1+a*b*gamma2+a*b*gamma3+b^2*gamma4]         ;
    optimize(Model,u)
    
end