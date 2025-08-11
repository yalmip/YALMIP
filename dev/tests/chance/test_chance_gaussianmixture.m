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

function test_gaussian_twoterms_two_components_merged(testCase)

% Case with two elementwise gaussian mixtures each with two components,
% meaning the linear component will be expanded to 2*2 = 4 components
% Vector w merged from indivual variables and thus merged mixtures
yalmip('clear');
gamma = sdpvar(1);
u = sdpvar(1);
w1 = sdpvar(1,1);
w2 = sdpvar(1,1);
a = 0.25;
b = 0.75;
P1 = probability(w1+w2<=u);
Model = [uncertain(w1,'normal mixture',{[-5],[5]},{[1],[sqrt(3);3]},{a, b}),
         uncertain(w2,'normal mixture',{[-3],[2]},{[2],[sqrt(3);3]},[a, b]),
    P1 >= 1-gamma, gamma <= 0.04]
optimize(Model,u)
testCase.assertTrue(abs(value(u)-12.0844) <= 1e-3)


function test_gaussian_twoterms_two_components(testCase)

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
    W11 = -5 + 1.*randn(1,a*N);
    W12 = 5 + sqrt(3).*randn(1,b*N);
    W21 = -3 + 2.*randn(1,a*N);
    W22 = 2 + 3.*randn(1,b*N);
    W = [W11 W12;W21 W22];    
    W(1,:) = W(1,randperm(length(W)));
    W(2,:) = W(2,randperm(length(W)));
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


function test_gaussian_twoterms_three_components(testCase)

% Case with two elementwise gaussian mixtures each with three components,
% meaning the linear component will be expanded to 2*3 = 6 components
yalmip('clear');
gamma = sdpvar(1);
u = sdpvar(1);
w = sdpvar(2,1);
a = 0.1;
b = 0.2;
c = 1-a-b;
P1 = probability(w(1)-2*w(2)<=u);
Model = [uncertain(w,'normal mixture',{[-5;-3],[0;0],[5;2]},{[1;2],[1;1],[sqrt(3);3]},{a, b, c}),
    P1 >= 1-gamma, gamma <= 0.04]
optimize(Model,u)
testCase.assertTrue(abs(value(u)-12.787) <= 1e-3)

% Numerical verification
if 0    
    N = 1e6;
    W11 = -5 + 1.*randn(1,a*N);
    W12 = 0 + 1.*randn(1,b*N);
    W13 = 5 + sqrt(3).*randn(1,c*N);
    W21 = -3 + 2.*randn(1,a*N);
    W22 = 0 + 1.*randn(1,b*N);
    W23 = 2 + 3.*randn(1,c*N);
    W = [W11 W12 W13;W21 W22 W23];    
    W(1,:) = W(1,randperm(length(W)));
    W(2,:) = W(2,randperm(length(W)));
    Y = [1 -2]*W;
    % Check
    [value(gamma) 1-nnz(Y < value(u))/length(Y)]                
end

function test_gaussian_twoterms_three_components_differentmix(testCase)

% Case with two elementwise gaussian mixtures each with three components,
% meaning the linear component will be expanded to 2*3 = 6 components
yalmip('clear');
gamma = sdpvar(1);
u = sdpvar(1);
w = sdpvar(2,1);
a = 0.1;
b = 0.2;
c = 1-a-b;
d = 0.5;
e = 0.45;
f = 1-d-e;
P1 = probability(w(1)-2*w(2)<=u);
Model = [uncertain(w,'normal mixture',{[-5;-3],[0;0],[5;2]},{[1;2],[1;1],[sqrt(3);3]},[a b c;d e f]),
    P1 >= 1-gamma, gamma <= 0.04]
optimize(Model,u)
testCase.assertTrue(abs(value(u)-16.3058) <= 1e-3)

% Numerical verification
if 0    
    N = 1e6;
    W11 = -5 + 1.*randn(1,a*N);
    W12 = 0 + 1.*randn(1,b*N);
    W13 = 5 + sqrt(3).*randn(1,c*N);
    W21 = -3 + 2.*randn(1,d*N);
    W22 = 0 + 1.*randn(1,e*N);
    W23 = 2 + 3.*randn(1,ceil(f*N));
    W = [W11 W12 W13;W21 W22 W23];    
    W(1,:) = W(1,randperm(length(W)));
    W(2,:) = W(2,randperm(length(W)));
    Y = [1 -2]*W;
    % Check
    [value(gamma) 1-nnz(Y < value(u))/length(Y)]                
end


function test_gaussian_threeterms_three_components(testCase)

% Case with three elementwise gaussian mixtures each with three components,
% meaning the linear component will be expanded to 3*3*3 = 27 components
yalmip('clear');
gamma = sdpvar(1);
u = sdpvar(1);
w = sdpvar(3,1);
a = 0.1;
b = 0.2;
c = 1-a-b;
P1 = probability(w(1)-2*w(2)+w(3)<=u);
Model = [uncertain(w,'normal mixture',{[-5;-3;-1],[0;0;0],[5;2;1]},{[1;2;1],[1;1;1],[sqrt(3);3;1]},{a, b, c}),
    P1 >= 1-gamma, gamma <= 0.04]
optimize(Model,u)
testCase.assertTrue(abs(value(u)-13.55) <= 1e-3)

% Numerical verification
if 0    
    N = 1e6;
    W11 = -5 + 1.*randn(1,a*N);
    W12 = 0 + 1.*randn(1,b*N);
    W13 = 5 + sqrt(3).*randn(1,c*N);
        
    W21 = -3 + 2.*randn(1,a*N);
    W22 = 0 + 1.*randn(1,b*N);
    W23 = 2 + 3.*randn(1,c*N);
    
    W31 = -1 + 1.*randn(1,a*N);
    W32 = 0 + 1.*randn(1,b*N);
    W33 = 1 + 1.*randn(1,c*N);
    
    W = [W11 W12 W13;W21 W22 W23;W31 W32 W33];    
    W(1,:) = W(1,randperm(length(W)));
    W(2,:) = W(2,randperm(length(W)));
    W(3,:) = W(3,randperm(length(W)));
    Y = [1 -2 1]*W;
    % Check
    [value(gamma) 1-nnz(Y < value(u))/length(Y)]                
end



function test_gaussian_threeterms_three_components_missing(testCase)

% Case with three elementwise gaussian mixtures each with three components,
% meaning the linear component will be expanded to 3^3 = 27 components, but
% one is not used to should expand to 3^2 components
yalmip('clear');
gamma = sdpvar(1);
u = sdpvar(1);
w = sdpvar(3,1);
a = 0.1;
b = 0.2;
c = 1-a-b;
P1 = probability(w(1)+w(3)<=u);
Model = [uncertain(w,'normal mixture',{[-5;-3;-1],[0;0;0],[5;2;1]},{[1;2;1],[1;1;1],[sqrt(3);3;1]},{a, b, c}),
    P1 >= 1-gamma, gamma <= 0.04]
optimize(Model,u)
testCase.assertTrue(abs(value(u)-8.902) <= 1e-3)

% Numerical verification
if 0    
    N = 1e6;
    W11 = -5 + 1.*randn(1,a*N);
    W12 = 0 + 1.*randn(1,b*N);
    W13 = 5 + sqrt(3).*randn(1,c*N);
        
    W21 = -3 + 2.*randn(1,a*N);
    W22 = 0 + 1.*randn(1,b*N);
    W23 = 2 + 3.*randn(1,c*N);
    
    W31 = -1 + 1.*randn(1,a*N);
    W32 = 0 + 1.*randn(1,b*N);
    W33 = 1 + 1.*randn(1,c*N);
    
    W = [W11 W12 W13;W21 W22 W23;W31 W32 W33];    
    W(1,:) = W(1,randperm(length(W)));
    W(2,:) = W(2,randperm(length(W)));
    W(3,:) = W(3,randperm(length(W)));
    Y = [1 0 1]*W;
    % Check
    [value(gamma) 1-nnz(Y < value(u))/length(Y)]                
end




function test_gaussian_decision_in_center(testCase)

yalmip('clear');
gamma = sdpvar(1);
u = sdpvar(1);
sdpvar m1 m2 m3 m4
w = sdpvar(2,1);
a = 0.25;
b = 0.75;

P1 = probability(w(1)+w(2)<=u)
P2 = probability(-u <= w(1)+w(2));
Model = [uncertain(w,'normal mixture',{[m1;m2],[m3;m4]},{[1;1],[1;1]},{a, b}),
    P1 >= 1-gamma, P2 >= 1-gamma, gamma <= 0.01]
optimize(Model,u)
testCase.assertTrue(norm(value([m1 m2 m3 m4])) <= 1e-3)
