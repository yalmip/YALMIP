function tests = test_chance_1
tests = functiontests(localfunctions);

function test4(testCase)
yalmip('clear')
w = sdpvar(2,1);
wmean = [3;4];
a = [1;3];
sdpvar s1 s2
Model = [uncertain(w,'exponential',wmean), probability(a'*w >= s1) >= .5,
                                           probability(a'*w <= s2) >= .5];
optimize(Model,abs(s1-s2))                   
testCase.assertTrue(abs(abs(value(s1-s2))-17) <= 6)

function test11(testCase)   
yalmip('clear')
x = sdpvar(4,1);
w1 = sdpvar(4,1);
w2 = sdpvar(4,1);
sdpvar s
Model = [uncertain(w1,'exponential',([1;2;3;4])),
         uncertain(w2,'exponential',([1;2;3;4])),
         probability(x'*(w1-w2) >= -1) >= s, sum(x)==1, x>=0];
optimize(Model,-s,sdpsettings('debug',1,'solver',''))

function test12(testCase)  
% Optimal sensor fusion of unkown with specified normal
yalmip('clear')
x = sdpvar(4,1);
w = sdpvar(4,1);
sdpvar s
Model = [uncertain(w,'moment',[0;0;0;0],diag([1;2;3;4])),
         probability(x'*w >= -1) >= s, sum(x)==1, x>=0];
optimize(Model,-s,sdpsettings('debug',1,'solver',''));
testCase.assertTrue(abs(value(s)-0.6756) <= 1e-3)

function test14(testCase)  
load portfoliodata
yalmip('clear')
w = sdpvar(N,1);
u = sdpvar(N,1);
alpha = sdpvar(1)
Markowitz = [0 <= w, sum(w) == 1];
VaRModel = [uncertain(u,'data',S), probability(u'*w >= -alpha) >= .95];
optimize([Markowitz, VaRModel], alpha)
testCase.assertTrue(abs(value(alpha)-0.0272) <= 1e-3)