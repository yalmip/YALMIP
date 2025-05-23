function tests = test_chance_1
tests = functiontests(localfunctions);

function test_gaussian_free_mean(testCase)

yalmip('clear')
sdpvar w wmean
Model = [uncertain(w,'normal', wmean, 1), 
         probability(w >= 0) >= 0.95];     
optimize(Model, wmean)
testCase.assertTrue(abs(value(wmean)-1.6449) <= 1e-3)

function test_gaussian_representmedian(testCase)
yalmip('clear')
w = sdpvar(2,1);
wmean = [3;4];
a = [1;3];
sdpvar s
Model = [uncertain(w,'normal',wmean,1), probability(a'*w >= s) >= .5,
                                        probability(a'*w <= s) >= .5];
optimize(Model)
testCase.assertTrue(abs(value(s)-15) <= 1e-3)
% nnz(a'*(wmean + randn(2,1e6))>=15)/1e6
% median(a'*(wmean + randn(2,1e6)))

function test_gaussian_worstcase_stddev(testCase)
yalmip('clear')
w=sdpvar(1,1);
sdpvar t

sdpvar R
Model = [probability(w >= t) >= 0.95,uncertain(w,'normal',0,1+R)];
Model = [Model, uncertain(R),0<=R<=1];
optimize(Model,-t)
testCase.assertTrue(abs(value(t)--3.2897) <= 1e-3)

function test_gaussian_worstcasemean(testCase)
yalmip('clear')
w = sdpvar(2,1);
wmean = sdpvar(2,1);
a = [1;3];
sdpvar s
Model = [uncertain(w,'normal',wmean,1),
         uncertain([2.9;3.9] <= wmean <= [3.1;4.1]), 
         probability(a'*w >= s) >= .5]
optimize(Model,-s)
testCase.assertTrue(abs(value(s)-14.6) <= 1e-3)

% Corresponding model
wmean = [2.9;3.9];
a = [1;3];
sdpvar s
Model = [uncertain(w,'normal',wmean,1),        
         probability(a'*w >= s) >= .5]
optimize(Model,-s)
testCase.assertTrue(abs(value(s)-14.6) <= 1e-3)

function test_gaussian_fusion(testCase)                                    
% Optimal sensor fusion of gaussians propotional to inverse variance

yalmip('clear')
x = sdpvar(4,1);
w = sdpvar(4,1);
sdpvar s
Model = [uncertain(w,'mvnrnd',[0;0;0;0],diag([1;2;3;4])),
         probability(x'*w >= -1) >= s, sum(x)==1, x>=0];
optimize(Model,-s)
testCase.assertTrue(all(abs(value(x)-[.48;.24;.16;.12]) <= 1e-3))

yalmip('clear')
x = sdpvar(4,1);
w = sdpvar(4,1);
sdpvar s
Model = [uncertain(w,'normal',[0;0;0;0],[1;2;3;4].^.5),
         probability(x'*w >= -1) >= s, sum(x)==1, x>=0];
optimize(Model,-s)
testCase.assertTrue(all(abs(value(x)-[.48;.24;.16;.12]) <= 1e-3))

x = sdpvar(4,1);
w = sdpvar(2,1);
q = sdpvar(2,1);
sdpvar s m
Model = [uncertain(w(1),'normal',0,1),        
         uncertain(w(2),'normal',0,sqrt(2)),                                             
         uncertain(q   ,'mvnrnd',[0;0],[3;4]),
            probability(x'*[w;q] >= -1) >= s, sum(x)==1, x>=0];
optimize(Model,-s)
testCase.assertTrue(all(abs(value(x)-[.48;.24;.16;.12]) <= 1e-3))

function test_gaussian_fusion_gaussian_mean(testCase)   
% Optimal sensor fusion of normal, merge distributions automatically
% Recursive as mean is normal, makes sense? Anyway, runs...
yalmip('clear')
x = sdpvar(4,1);
w = sdpvar(2,1);
q = sdpvar(2,1);
sdpvar s m1 m2 m3 m4
Model = [uncertain(w(1),'normal',m1,1),
         uncertain(w(2),'normal',m2,sqrt(2)),         
         uncertain(q ,'normal',[m3;m4],[3;4].^.5),         
         uncertain(m1,'normal',0,.5.^.5),
         uncertain(m2,'normal',0,.6.^.5),
         uncertain(m3,'normal',0,.7.^.5),
         uncertain(m4,'normal',0,.8.^.5),
         probability(x'*[w;q] >= -1) >= s, sum(x)==1, x>=0];
optimize(Model,-s)
testCase.assertTrue(abs(value(s)-.8169) <= 1e-3)


function test_gaussian_value_at_risk(testCase)   
load portfoliodata
yalmip('clear')
w = sdpvar(N,1);
u = sdpvar(N,1);
alpha = sdpvar(1)
Markowitz = [0 <= w, sum(w) == 1];
VaRModel = [uncertain(u,'mvnrnd',r,R), probability(w'*u >= -alpha) >= .95];
optimize([Markowitz, VaRModel], alpha)
testCase.assertTrue(abs(value(alpha)-0.0189) <= 1e-3)

