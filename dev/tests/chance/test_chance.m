function tests = test_chance
tests = functiontests(localfunctions);

function test_simple_normal(testCase)

yalmip('clear')
a=sdpvar(1,1);
sdpvar t

Model = [probability(a >= t) >= 0.5,uncertain(a,'normal',[0],eye(1))];
optimize(Model,-t)
testCase.assertTrue(abs(value(t)) <= 1e-4)

Model = [probability(a >= t) >= 0.5,uncertain(a,'normalf',[0],eye(1))];
optimize(Model,-t)
testCase.assertTrue(abs(value(t)) <= 1e-4)

Model = [probability(a >= t) >= 0.95,uncertain(a,'normal',[0],4)];
optimize(Model,-t)
testCase.assertTrue(abs(value(t)--3.2897) <= 1e-4)

Model = [probability(a >= t) >= 0.95,uncertain(a,'normalf',[0],2)];
optimize(Model,-t)
testCase.assertTrue(abs(value(t)--3.2897) <= 1e-4)


function test_robust_normal(testCase)
yalmip('clear')
a=sdpvar(1,1);
sdpvar t

% worst case should be R = 1
sdpvar R
Model = [probability(a >= t) >= 0.95,uncertain(a,'normalf',[0],1+R)];
Model = [Model, uncertain(R),0<=R<=1];
optimize(Model,-t)
testCase.assertTrue(abs(value(t)--3.2897) <= 1e-4)



function test_complex_normal(testCase)

yalmip('clear')
sdpvar a b t
sdpvar t

Model = [probability(a >= t) >= 0.5,
         probability(b >= t) >= 0.5, 
         uncertain(a,'normal',[0],eye(1)),
         uncertain(b,'normal',[0],eye(1))];
optimize(Model,-t)
testCase.assertTrue(abs(value(t)-0) <= 1e-4)

sdpvar w wmean
Model = [uncertain(w,'normal', wmean, 1), 
         probability(w >= 0) >= 0.9];     
optimize(Model, wmean)
testCase.assertTrue(abs(value(wmean)-1.28155) <= 1e-4)


%yalmip('clear')
%sdpvar a b t m
%Model = [probability(a >= t) >= 0.5,
%         probability(b >= t) >= 0.5, 
%         uncertain(a,'normal',[0],eye(1)),
%         uncertain(b,'normal',m,eye(1)),
%         uncertain(m,'normal',0,eye(1)), ];
%P=optimizer(Model,-t,sdpsettings('solver','cplex'),m,t)
%Q=sample(P,20)

w = sdpvar(2,1);
wmean = [3;4];
a = [1;3];
sdpvar s
Model = [uncertain(w,'normal',wmean,1), probability(a'*w >= s) >= .5,
                                        probability(a'*w <= s) >= .5];
optimize(Model)
testCase.assertTrue(abs(value(s)-15) <= 1e-4)

%mbg_asserttolequal(double(s), 15, 1e-3);

sdpvar s1 s2
Model = [uncertain(w,'exponential',wmean), probability(a'*w >= s1) >= .5,
                                           probability(a'*w <= s2) >= .5];
optimize(Model,abs(s1-s2))                                   
mbg_asserttolequal(value(abs(s1-s2)), 17, 6);

w = sdpvar(2,1);
wmean = [3;4];
a = [1;3];
sdpvar s

Model = [uncertain(w,'normal',wmean,1), probability(a'*w >= s) >= .9]
optimize(Model,-s)
mbg_asserttolequal(value(s), 10.94, 1e-2);


w = sdpvar(2,1);
wmean = sdpvar(2,1);
a = [1;3];
sdpvar s

Model = [uncertain(w,'normal',wmean,1),
         uncertain([2.9;3.9] <= wmean <= [3.1;4.1]), 
         probability(a'*w >= s) >= .5]
optimize(Model,-s)
mbg_asserttolequal(value(s), 14.6, 1e-2);

                                    
                
% Optimal sensor fusion of normal
clear all
x = sdpvar(4,1);
w = sdpvar(4,1);
sdpvar s
Model = [uncertain(w,'normalm',[0;0;0;0],diag([1;2;3;4])),
         probability(x'*w >= -1) >= s, sum(x)==1, x>=0];
optimize(Model,-s,sdpsettings('debug',1,'solver',''))
mbg_asserttolequal(value(x), [.48;.24;.16;.12], 1e-2);

% Optimal sensor fusion of normal
clear all
x = sdpvar(4,1);
w = sdpvar(4,1);
sdpvar s
Model = [uncertain(w,'normal',[0;0;0;0],([1;2;3;4])),
         probability(x'*w >= -1) >= s, sum(x)==1, x>=0];
optimize(Model,-s,sdpsettings('debug',1,'solver',''))
mbg_asserttolequal(value(x), [.48;.24;.16;.12], 1e-2);


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
mbg_asserttolequal(value(x), [.48;.24;.16;.12], 1e-2);


% Optimal sensor fusion of normal, merge distributions automatically
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
mbg_asserttolequal(value(s), .8169, 1e-2);

clear all
x = sdpvar(4,1);
w1 = sdpvar(4,1);
w2 = sdpvar(4,1);
sdpvar s
Model = [uncertain(w1,'exponential',([1;2;3;4])),
         uncertain(w2,'exponential',([1;2;3;4])),
         probability(x'*(w1-w2) >= -1) >= s, sum(x)==1, x>=0];
optimize(Model,-s,sdpsettings('debug',1,'solver',''))

% Optimal sensor fusion of unkown with specified normal
clear all
x = sdpvar(4,1);
w = sdpvar(4,1);
sdpvar s
Model = [uncertain(w,'moment',[0;0;0;0],diag([1;2;3;4])),
         probability(x'*w >= -1) >= s, sum(x)==1, x>=0];
optimize(Model,-s,sdpsettings('debug',1,'solver',''));
mbg_asserttolequal(value(s), 0.6756, 1e-2);

load portfoliodata
w = sdpvar(N,1);
u = sdpvar(N,1);
alpha = sdpvar(1)
Markowitz = [0 <= w, sum(w) == 1];
VaRModel = [uncertain(u,'normalm',r,R), probability(w'*u >= -alpha) >= .95];
optimize([Markowitz, VaRModel], alpha)
LowestVaR = value(alpha);
mbg_asserttolequal(value(alpha), .0189, 1e-2);

Markowitz = [0 <= w, sum(w) == 1];
VaRModel = [uncertain(u,'data',S), probability(u'*w >= -alpha) >= .95];
optimize([Markowitz, VaRModel], alpha)
mbg_asserttolequal(value(alpha), 0.0272, 1e-2);