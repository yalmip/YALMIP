function tests = test_robust_5
tests = functiontests(localfunctions);

function test1(testCase)
% Problem with two separable polytope uncertainties
randn('seed',0);
A1 = randn(8,2);
A2 = randn(8,2);
x0 = randn(2,1);
b1 = A1*x0+1;
b2 = A2*x0+1;

w2 = sdpvar(2,1);
w1 = sdpvar(2,1);
x = sdpvar(1,1);

sol1 = optimize([uncertain([w1;w2]), x+sum(w1)+sum(w2)<=0,A1*w1 <= b1, A2*w2 <= b2],-x,sdpsettings('verbose',0));
obj1 = value(-x);
sol2 = optimize([A1*w1 <= b1,A2*w2 <= b2],-sum(w1)-sum(w2),sdpsettings('verbose',0));
obj2 = value(-sum(w1)-sum(w2));

testCase.assertTrue(sol1.problem == 0);
testCase.assertTrue(abs(obj1+obj2) <= 1e-4)

function test2(testCase)
% 2 balls
x = sdpvar(1);
w = sdpvar(3,1); 
q = sdpvar(2,1);
sol = optimize([x+sum(w)+sum(q)<=1,-1<=w<=3,uncertain(w),uncertain(q),q'*q<=2],-x,sdpsettings('verbose',0));
testCase.assertTrue(sol.problem == 0);
testCase.assertTrue(abs(value(x)--10) <= 1e-5)

function test3(testCase)
% 2 polytopic
x = sdpvar(1);
w = sdpvar(3,1); 
q = sdpvar(2,1);
sol = optimize([x+sum(w)+sum(q)<=1,-1<=w<=3,uncertain(w),uncertain(q),q'*q<=2,q>=0,w>=0],-x,sdpsettings('verbose',0));
testCase.assertTrue(sol.problem == 0);
testCase.assertTrue(abs(value(x)--10) <= 1e-5)

function test4(testCase)
% Messy with uncertain auxilliary variables
yalmip('clear')
x = sdpvar(1);
w1 = sdpvar(1);
w2 = sdpvar(1);
t = sdpvar(1);
A = [1 1;-1 0;0 -2;-1 -1;0 -1/2];
b = [.1;.1;.1;.1;.1];
c = randn(2,1);
ops = sdpsettings('verbose',0,'robust.auxreduce','enumeration');
optimize([x+max(x+c'*[w1;w2],norm([x*w1 + A*[w1;w2]],2))+norm([x - A*[w1;w2]],1)<=t,uncertain([w1;w2]),A*[w1;w2] <= b],t,ops)
o1 = value(t);

v = [-.1 -.05 .15 -.1;
     0 -0.05 -0.05 0.2]';
C = [];
for i = 1:size(v,1)
C = [C,x+max(x+c'*v(i,:)',norm(x*v(i,1) + A*v(i,:)',2))+norm(x - A*v(i,:)',1)<=t];
end
optimize(C,t,sdpsettings('verbose',0));
o2 = value(t);
testCase.assertTrue(abs(o1-o2) <= 1e-5)


function test5(testCase)
% Messy with uncertain auxilliary variables
yalmip('clear')
x = sdpvar(1);
w1 = sdpvar(1);
w2 = sdpvar(1);
t = sdpvar(1);
A = [1 1;-1 0;0 -2;-1 -1;0 -1/2];
b = [.1;.1;.1;.1;.1];
c = randn(2,1);
ops = sdpsettings('verbose',1,'robust.auxreduce','enumeration');
optimize([x+max(x+c'*[w1;w2],norm([x*w1 + A*[w1;w2]],1))+norm([x - A*[w1;w2]],1)<=t,uncertain([w1;w2]),A*[w1;w2] <= b],t,ops);
o1 = value(t);

v = [-.1 -.05 .15 -.1;
     0 -0.05 -0.05 0.2]';
C = [];
for i = 1:size(v,1)
C = [C,x+max(x+c'*v(i,:)',norm(x*v(i,1) + A*v(i,:)',1))+norm(x - A*v(i,:)',1)<=t];
end
optimize(C,t,sdpsettings('verbose',0));
o2 = value(t);
testCase.assertTrue(abs(o1-o2) <= 1e-5)

function test6(testCase)
% Should yield 27
yalmip('clear')
x = sdpvar(1);
w = sdpvar(2,1);
C = [x+norm(w+1,1) <= 30];
W = [uncertain(w),norm(w+.5,1) + norm(w,1)<=3];
ops = sdpsettings('verbose',0,'robust.auxreduce','enumeration');
optimize([C,W],-x,ops);
testCase.assertTrue(abs(value(x) - 27) <= 1e-5)

function test7(testCase)
% Simple explicit maximization
yalmip('clear')
x = sdpvar(1);
w = sdpvar(2,1);
C = [x+w(1)+w(2) <= 1];
W = [uncertain(w),norm(w,1)<=2];
optimize([C,W],-x,sdpsettings('verbose',0));
testCase.assertTrue(abs(value(x)--1) <= 1e-5)

function test8(testCase)
yalmip('clear')
x = sdpvar(1);
w = sdpvar(2,1);
C = [x+w(1)+w(2) <= 1];
W = [uncertain(w),norm(w,2)<=2];
optimize([C,W],-x,sdpsettings('verbose',0));
testCase.assertTrue(abs(value(x)--1.8284) <= 1e-2)

function test9(testCase)
yalmip('clear')
x = sdpvar(1);
w = sdpvar(2,1);
C = [x+w(1)+w(2) <= 1];
W = [uncertain(w),norm(w,inf)<=2];
optimize([C,W],-x,sdpsettings('verbose',0));
testCase.assertTrue(abs(value(x)--3)<=1e-2)

function test10(testCase)
% Missing term in explicit expression
yalmip('clear')
x = sdpvar(1);
w = sdpvar(2,1);
C = [x+2*w(2) <= 1];
W = [uncertain(w),norm(w,1)<=2];
optimize([C,W],-x,sdpsettings('verbose',0));
testCase.assertTrue(abs(value(x)--3)<=1e-5)

function test11(testCase)
yalmip('clear')
x = sdpvar(1);
w = sdpvar(2,1);
C = [x+2*w(1) <= 1];
W = [uncertain(w),norm(w,2)<=2];
optimize([C,W],-x,sdpsettings('verbose',0));
testCase.assertTrue(abs(value(x)--3)<=1e-5)

function test12(testCase)
yalmip('clear')
x = sdpvar(1);
w = sdpvar(2,1);
C = [x+2*w(1) <= 1];
W = [uncertain(w),norm(w,inf)<=2];
optimize([C,W],-x,sdpsettings('verbose',0));
testCase.assertTrue(abs(value(x)--3)<=1e-5)

function test13(testCase)
% Enumeration of 6 vertices
yalmip('clear')
x = sdpvar(1);
w = sdpvar(2,1);
C = [x+norm(w,2) <= 1];
W = [uncertain(w),norm(w,1)<=2];
optimize([C,W],-x,sdpsettings('verbose',1,'solver','sedumi'))
value(x)
testCase.assertTrue((abs(value(x)--1)<=1e-1) || (sol.problem == 4))

function test14(testCase)
% Projection + explicit maximization
yalmip('clear')
x = sdpvar(1);
w = sdpvar(2,1);
C = [x+norm(w,1) <= 1];
W = [uncertain(w),norm(w,2)<=2];
optimize([C,W],-x,sdpsettings('verbose',0,'robust.auxreduce','projection'));
testCase.assertTrue(abs(value(x)--1.8284) <= 1e-3)

function test15(testCase)
yalmip('clear')
x = sdpvar(1);
w = sdpvar(2,1);
C = [x+norm(w,2) <= 1];
W = [uncertain(w),norm(w,2)<=2];
optimize([C,W],-x,sdpsettings('verbose',1,'robust.auxreduce','projection','sedumi.free',0));
testCase.assertTrue(abs(value(x)--1) <= 1e-1)

function test16(testCase)
yalmip('clear')
x = sdpvar(1);
w = sdpvar(2,1);
C = [x+norm(w,2)+norm(w,1) <= 1];
W = [uncertain(w),norm(w,2)+norm(w,1)<=2];
optimize([C,W],-x,sdpsettings('verbose',0,'robust.auxreduce','projection','sedumi.free',0));
testCase.assertTrue(abs(value(x)--1) <= 1e-1)

function test17(testCase)
yalmip('clear')
x = sdpvar(1);
w = sdpvar(2,1);
C = [x+norm(w,1)+norm(w,2) <= 1];
W = [uncertain(w),norm(w,2)<=2];
optimize([C,W],-x,sdpsettings('verbose',0,'robust.auxreduce','projection'));
testCase.assertTrue(abs(value(x)--3.8284) <= 1e-1)

function test18(testCase)
x = sdpvar(1);
w = sdpvar(1);
C = [1<=x<=4];
W = [uncertain(w),-10<=w<=2];
sol = optimize([C,W],x*w^2,sdpsettings('verbose',0));
testCase.assertTrue(sol.problem==1 | sol.problem == 12)

C = [1<=x<=4];
W = [uncertain(w),-10<=w<=2];
sol = optimize([C,W],x*w,sdpsettings('verbose',0));
testCase.assertTrue(abs(value(x)-1)<=1e-1)