function test_robust_5


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

sol1 = solvesdp([uncertain([w1;w2]), x+sum(w1)+sum(w2)<=0,A1*w1 <= b1, A2*w2 <= b2],-x,sdpsettings('debug',1))
obj1 = double(-x);
sol2 = solvesdp([A1*w1 <= b1,A2*w2 <= b2],-sum(w1)-sum(w2))
obj2 = double(-sum(w1)-sum(w2));

mbg_asserttolequal(sol1.problem, 0, 1e-4);
mbg_asserttolequal(obj1+obj2,0, 1e-4);

% 2 balls
x = sdpvar(1);
w = sdpvar(3,1); 
q = sdpvar(2,1);
sol = solvesdp([x+sum(w)+sum(q)<=1,-1<=w<=3,uncertain(w),uncertain(q),q'*q<=2],-x)
mbg_asserttolequal(sol.problem, 0, 1e-5);
mbg_asserttolequal(double(x),-10,1e-5);

% 2 polytopic
x = sdpvar(1);
w = sdpvar(3,1); 
q = sdpvar(2,1);
sol = solvesdp([x+sum(w)+sum(q)<=1,-1<=w<=3,uncertain(w),uncertain(q),q'*q<=2,q>=0,w>=0],-x)
mbg_asserttolequal(sol.problem, 0, 1e-5);
mbg_asserttolequal(double(x),-10, 1e-5);

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
solvesdp([x+max(x+c'*[w1;w2],norm([x*w1 + A*[w1;w2]],2))+norm([x - A*[w1;w2]],1)<=t,uncertain([w1;w2]),A*[w1;w2] <= b],t,ops)
o1 = double(t);

v = extreme(polytope(A,b));
C = [];
for i = 1:size(v,1)
C = [C,x+max(x+c'*v(i,:)',norm(x*v(i,1) + A*v(i,:)',2))+norm(x - A*v(i,:)',1)<=t];
end
solvesdp(C,t)
o2 = double(t);
mbg_asserttolequal(o1-o2,0, 1e-5);



% Messy with uncertain auxilliary variables
yalmip('clear')
x = sdpvar(1);
w1 = sdpvar(1);
w2 = sdpvar(1);
t = sdpvar(1);
A = [1 1;-1 0;0 -2;-1 -1;0 -1/2];
b = [.1;.1;.1;.1;.1];
c = randn(2,1);
solvesdp([x+max(x+c'*[w1;w2],norm([x*w1 + A*[w1;w2]],1))+norm([x - A*[w1;w2]],1)<=t,uncertain([w1;w2]),A*[w1;w2] <= b],t,ops);
o1 = double(t);

v = extreme(polytope(A,b));
C = [];
for i = 1:size(v,1)
C = [C,x+max(x+c'*v(i,:)',norm(x*v(i,1) + A*v(i,:)',1))+norm(x - A*v(i,:)',1)<=t];
end
solvesdp(C,t)
o2 = double(t);
mbg_asserttolequal(o1-o2,0, 1e-5);


% Should yield 27
yalmip('clear')
x = sdpvar(1);
w = sdpvar(2,1);
C = [x+norm(w+1,1) <= 30]
W = [uncertain(w),norm(w+.5,1) + norm(w,1)<=3];
solvesdp([C,W],-x,ops);
mbg_asserttolequal(double(x),27,1e-5)


% Simple explicit maximization
yalmip('clear')
x = sdpvar(1);
w = sdpvar(2,1);
C = [x+w(1)+w(2) <= 1]
W = [uncertain(w),norm(w,1)<=2];
solvesdp([C,W],-x)
mbg_asserttolequal(double(x),-1,1e-5)

yalmip('clear')
x = sdpvar(1);
w = sdpvar(2,1);
C = [x+w(1)+w(2) <= 1]
W = [uncertain(w),norm(w,2)<=2];
solvesdp([C,W],-x)
mbg_asserttolequal(double(x),-1.8284,1e-2)

yalmip('clear')
x = sdpvar(1);
w = sdpvar(2,1);
C = [x+w(1)+w(2) <= 1]
W = [uncertain(w),norm(w,inf)<=2];
solvesdp([C,W],-x)
mbg_asserttolequal(double(x),-3,1e-2)

% Missing term in explicit expression
yalmip('clear')
x = sdpvar(1);
w = sdpvar(2,1);
C = [x+2*w(2) <= 1]
W = [uncertain(w),norm(w,1)<=2];
solvesdp([C,W],-x)
mbg_asserttolequal(double(x),-3,1e-5)
yalmip('clear')
x = sdpvar(1);
w = sdpvar(2,1);
C = [x+2*w(1) <= 1]
W = [uncertain(w),norm(w,2)<=2];
solvesdp([C,W],-x)
mbg_asserttolequal(double(x),-3,1e-5)
yalmip('clear')
x = sdpvar(1);
w = sdpvar(2,1);
C = [x+2*w(1) <= 1]
W = [uncertain(w),norm(w,inf)<=2];
solvesdp([C,W],-x)
mbg_asserttolequal(double(x),-3,1e-5)


% Enumeration of 6 vertices
yalmip('clear')
x = sdpvar(1);
w = sdpvar(2,1);
C = [x+norm(w,2) <= 1]
W = [uncertain(w),norm(w,1)<=2];
solvesdp([C,W],-x)
mbg_asserttolequal(double(x),-1,1e-1)


% Projection + explicit maximization
yalmip('clear')
x = sdpvar(1);
w = sdpvar(2,1);
C = [x+norm(w,1) <= 1]
W = [uncertain(w),norm(w,2)<=2];
solvesdp([C,W],-x,sdpsettings('robust.auxreduce','projection'))
mbg_asserttolequal(double(x),-1.8284,1e-3)


yalmip('clear')
x = sdpvar(1);
w = sdpvar(2,1);
C = [x+norm(w,2) <= 1]
W = [uncertain(w),norm(w,2)<=2];
solvesdp([C,W],-x,sdpsettings('robust.auxreduce','projection'))
mbg_asserttolequal(double(x),-1,1e-1)

yalmip('clear')
x = sdpvar(1);
w = sdpvar(2,1);
C = [x+norm(w,2)+norm(w,1) <= 1]
W = [uncertain(w),norm(w,2)+norm(w,1)<=2];
solvesdp([C,W],-x,sdpsettings('robust.auxreduce','projection'))
mbg_asserttolequal(double(x),-1,1e-1)

yalmip('clear')
x = sdpvar(1);
w = sdpvar(2,1);
C = [x+norm(w,1)+norm(w,2) <= 1]
W = [uncertain(w),norm(w,2)<=2];
solvesdp([C,W],-x,sdpsettings('robust.auxreduce','projection'))
mbg_asserttolequal(double(x),-3.8284,1e-1)

x = sdpvar(1);
w = sdpvar(1);
C = [1<=x<=4];
W = [uncertain(w),-10<=w<=2];
sol = solvesdp([C,W],x*w^2)
mbg_asserttrue(sol.problem==1 | sol.problem == 12)

C = [1<=x<=4];
W = [uncertain(w),-10<=w<=2];
sol = solvesdp([C,W],x*w)
mbg_asserttolequal(double(x),1,1e-1)


