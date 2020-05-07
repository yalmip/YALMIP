function tests = test_robust_5
tests = functiontests(localfunctions);

function test1(dummy)
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

sol1 = optimize([uncertain([w1;w2]), x+sum(w1)+sum(w2)<=0,A1*w1 <= b1, A2*w2 <= b2],-x,sdpsettings('debug',1))
obj1 = value(-x);
sol2 = optimize([A1*w1 <= b1,A2*w2 <= b2],-sum(w1)-sum(w2))
obj2 = value(-sum(w1)-sum(w2));

assert(sol1.problem == 0);
assert(abs(obj1+obj2) <= 1e-4)

function test2(dummy)
% 2 balls
x = sdpvar(1);
w = sdpvar(3,1); 
q = sdpvar(2,1);
sol = optimize([x+sum(w)+sum(q)<=1,-1<=w<=3,uncertain(w),uncertain(q),q'*q<=2],-x)
assert(sol.problem == 0);
assert(abs(value(x)--10) <= 1e-5)

function test3(dummy)
% 2 polytopic
x = sdpvar(1);
w = sdpvar(3,1); 
q = sdpvar(2,1);
sol = optimize([x+sum(w)+sum(q)<=1,-1<=w<=3,uncertain(w),uncertain(q),q'*q<=2,q>=0,w>=0],-x)
assert(sol.problem == 0);
assert(abs(value(x)--10) <= 1e-5)

function test4(dummy)
% Messy with uncertain auxilliary variables
yalmip('clear')
x = sdpvar(1);
w1 = sdpvar(1);
w2 = sdpvar(1);
t = sdpvar(1);
A = [1 1;-1 0;0 -2;-1 -1;0 -1/2];
b = [.1;.1;.1;.1;.1];
c = randn(2,1);
ops = sdpsettings('robust.auxreduce','enumeration');
optimize([x+max(x+c'*[w1;w2],norm([x*w1 + A*[w1;w2]],2))+norm([x - A*[w1;w2]],1)<=t,uncertain([w1;w2]),A*[w1;w2] <= b],t,ops)
o1 = value(t);

v = extreme(polytope(A,b));
C = [];
for i = 1:size(v,1)
C = [C,x+max(x+c'*v(i,:)',norm(x*v(i,1) + A*v(i,:)',2))+norm(x - A*v(i,:)',1)<=t];
end
optimize(C,t)
o2 = value(t);
assert(abs(o1-o2) <= 1e-5)


function test5(dummy)
% Messy with uncertain auxilliary variables
yalmip('clear')
x = sdpvar(1);
w1 = sdpvar(1);
w2 = sdpvar(1);
t = sdpvar(1);
A = [1 1;-1 0;0 -2;-1 -1;0 -1/2];
b = [.1;.1;.1;.1;.1];
c = randn(2,1);
ops = sdpsettings('robust.auxreduce','enumeration');
optimize([x+max(x+c'*[w1;w2],norm([x*w1 + A*[w1;w2]],1))+norm([x - A*[w1;w2]],1)<=t,uncertain([w1;w2]),A*[w1;w2] <= b],t,ops);
o1 = value(t);

v = extreme(polytope(A,b));
C = [];
for i = 1:size(v,1)
C = [C,x+max(x+c'*v(i,:)',norm(x*v(i,1) + A*v(i,:)',1))+norm(x - A*v(i,:)',1)<=t];
end
optimize(C,t)
o2 = value(t);
assert(abs(o1-o2) <= 1e-5)

function test6(dummy)
% Should yield 27
yalmip('clear')
x = sdpvar(1);
w = sdpvar(2,1);
C = [x+norm(w+1,1) <= 30]
W = [uncertain(w),norm(w+.5,1) + norm(w,1)<=3];
ops = sdpsettings('robust.auxreduce','enumeration');
optimize([C,W],-x,ops);
assert(abs(value(x) - 27) <= 1e-5)

function test7(dummy)
% Simple explicit maximization
yalmip('clear')
x = sdpvar(1);
w = sdpvar(2,1);
C = [x+w(1)+w(2) <= 1]
W = [uncertain(w),norm(w,1)<=2];
optimize([C,W],-x)
assert(abs(value(x)--1) <= 1e-5)

function test8(dummy)
yalmip('clear')
x = sdpvar(1);
w = sdpvar(2,1);
C = [x+w(1)+w(2) <= 1]
W = [uncertain(w),norm(w,2)<=2];
optimize([C,W],-x)
assert(abs(value(x)--1.8284) <= 1e-2)

function test9(dummy)
yalmip('clear')
x = sdpvar(1);
w = sdpvar(2,1);
C = [x+w(1)+w(2) <= 1]
W = [uncertain(w),norm(w,inf)<=2];
optimize([C,W],-x)
assert(abs(value(x)--3)<=1e-2)

function test10(dummy)
% Missing term in explicit expression
yalmip('clear')
x = sdpvar(1);
w = sdpvar(2,1);
C = [x+2*w(2) <= 1]
W = [uncertain(w),norm(w,1)<=2];
optimize([C,W],-x)
assert(abs(value(x)--3)<=1e-5)

function test11(dummy)
yalmip('clear')
x = sdpvar(1);
w = sdpvar(2,1);
C = [x+2*w(1) <= 1]
W = [uncertain(w),norm(w,2)<=2];
optimize([C,W],-x)
assert(abs(value(x)--3)<=1e-5)

function test12(dummy)
yalmip('clear')
x = sdpvar(1);
w = sdpvar(2,1);
C = [x+2*w(1) <= 1]
W = [uncertain(w),norm(w,inf)<=2];
optimize([C,W],-x)
assert(abs(value(x)--3)<=1e-5)

function test13(dummy)
% Enumeration of 6 vertices
yalmip('clear')
x = sdpvar(1);
w = sdpvar(2,1);
C = [x+norm(w,2) <= 1]
W = [uncertain(w),norm(w,1)<=2];
optimize([C,W],-x)
assert(abs(value(x)--1)<=1e-1)

function test14(dummy)
% Projection + explicit maximization
yalmip('clear')
x = sdpvar(1);
w = sdpvar(2,1);
C = [x+norm(w,1) <= 1]
W = [uncertain(w),norm(w,2)<=2];
optimize([C,W],-x,sdpsettings('robust.auxreduce','projection'))
assert(abs(value(x)--1.8284) <= 1e-3)

function test15(dummy)
yalmip('clear')
x = sdpvar(1);
w = sdpvar(2,1);
C = [x+norm(w,2) <= 1]
W = [uncertain(w),norm(w,2)<=2];
optimize([C,W],-x,sdpsettings('robust.auxreduce','projection'))
assert(abs(value(x)--1) <= 1e-1)

function test16(dummy)
yalmip('clear')
x = sdpvar(1);
w = sdpvar(2,1);
C = [x+norm(w,2)+norm(w,1) <= 1]
W = [uncertain(w),norm(w,2)+norm(w,1)<=2];
optimize([C,W],-x,sdpsettings('robust.auxreduce','projection'))
assert(abs(value(x)--1) <= 1e-1)

function test17(dummy)
yalmip('clear')
x = sdpvar(1);
w = sdpvar(2,1);
C = [x+norm(w,1)+norm(w,2) <= 1]
W = [uncertain(w),norm(w,2)<=2];
optimize([C,W],-x,sdpsettings('robust.auxreduce','projection'))
assert(abs(value(x)--3.8284) <= 1e-1)

function test18(dummy)
x = sdpvar(1);
w = sdpvar(1);
C = [1<=x<=4];
W = [uncertain(w),-10<=w<=2];
sol = optimize([C,W],x*w^2)
assert(sol.problem==1 | sol.problem == 12)

C = [1<=x<=4];
W = [uncertain(w),-10<=w<=2];
sol = optimize([C,W],x*w)
assert(abs(value(x)-1)<=1e-1)