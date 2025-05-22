function tests = test_chance_1
tests = functiontests(localfunctions);

function test1(testCase)
yalmip('clear')
sdpvar a b t
sdpvar t

Model = [probability(a >= t) >= 0.5,
         probability(b >= t) >= 0.5, 
         uncertain(a,'normal',[0],eye(1)),
         uncertain(b,'normal',[0],eye(1))];
optimize(Model,-t)
testCase.assertTrue(abs(value(t)) <= 1e-3)

function test2(testCase)
yalmip('clear')
sdpvar w wmean
Model = [uncertain(w,'normal', wmean, 1), 
         probability(w >= 0) >= 0.9];     
optimize(Model, wmean,sdpsettings('solver','','debug',1))
testCase.assertTrue(abs(value(wmean)-1.2816) <= 1e-3)


function test3(testCase)
yalmip('clear')
w = sdpvar(2,1);
wmean = [3;4];
a = [1;3];
sdpvar s
Model = [uncertain(w,'normal',wmean,1), probability(a'*w >= s) >= .5,
                                        probability(a'*w <= s) >= .5];
optimize(Model)
testCase.assertTrue(abs(value(s)-15) <= 1e-3)

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

function test5(testCase)
yalmip('clear')
w = sdpvar(2,1);
wmean = [3;4];
a = [1;3];
sdpvar s
Model = [uncertain(w,'normal',wmean,1), probability(a'*w >= s) >= .9]
optimize(Model,-s)
testCase.assertTrue(abs(value(s)-10.9474) <= 1e-3)

function test6(testCase)
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

function test7(testCase)                                    
% Optimal sensor fusion of normal
yalmip('clear')
x = sdpvar(4,1);
w = sdpvar(4,1);
sdpvar s
Model = [uncertain(w,'normalm',[0;0;0;0],diag([1;2;3;4])),
         probability(x'*w >= -1) >= s, sum(x)==1, x>=0];
optimize(Model,-s,sdpsettings('debug',1,'solver',''))
testCase.assertTrue(all(abs(value(x)-[.48;.24;.16;.12]) <= 1e-3))

function test8(testCase)                                    
% Optimal sensor fusion of normal
yalmip('clear')
x = sdpvar(4,1);
w = sdpvar(4,1);
sdpvar s
Model = [uncertain(w,'normal',[0;0;0;0],([1;2;3;4])),
         probability(x'*w >= -1) >= s, sum(x)==1, x>=0];
optimize(Model,-s,sdpsettings('debug',1,'solver',''))
testCase.assertTrue(all(abs(value(x)-[.48;.24;.16;.12]) <= 1e-3))

function test9(testCase)   
% Optimal sensor fusion of normal, merge distributions automatically
yalmip('clear')
x = sdpvar(4,1);
w = sdpvar(2,1);
q = sdpvar(2,1);
sdpvar s m
Model = [   uncertain(w(1),'normal',0,[1]),        
            uncertain(w(2),'normalf',0,[sqrt(2)]),                                             
            uncertain(q ,'normal',[0;0],[3;4]),
            probability(x'*[w;q] >= -1) >= s, sum(x)==1, x>=0];
optimize(Model,-s,sdpsettings('debug',1,'solver',''))
testCase.assertTrue(all(abs(value(x)-[.48;.24;.16;.12]) <= 1e-3))

function test10(testCase)   
% Optimal sensor fusion of normal, merge distributions automatically
% Recursive as mean is normal...
yalmip('clear')
x = sdpvar(4,1);
w = sdpvar(2,1);
q = sdpvar(2,1);
sdpvar s m1 m2 m3 m4
Model = [uncertain(w(1),'normal',m1,1),
         uncertain(w(2),'normal',m2,2),         
         uncertain(q ,'normal',[m3;m4],[3;4]),         
         uncertain(m1,'normal',0,.5),
         uncertain(m2,'normal',0,.6),
         uncertain(m3,'normal',0,.7),
         uncertain(m4,'normal',0,.8),
         probability(x'*[w;q] >= -1) >= s, sum(x)==1, x>=0];
optimize(Model,-s,sdpsettings('debug',1,'solver',''))
testCase.assertTrue(abs(value(s)-.8169) <= 1e-3)

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

function test13(testCase)  
load portfoliodata
yalmip('clear')
w = sdpvar(N,1);
u = sdpvar(N,1);
alpha = sdpvar(1)
Markowitz = [0 <= w, sum(w) == 1];
VaRModel = [uncertain(u,'normalm',r,R), probability(w'*u >= -alpha) >= .95];
optimize([Markowitz, VaRModel], alpha)
testCase.assertTrue(abs(value(alpha)-0.0189) <= 1e-3)

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