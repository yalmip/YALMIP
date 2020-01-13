function tests = test_operator_implies
tests = functiontests(localfunctions);

function test1(dummy)

% Binary variable implies LP constriants
sdpvar y u
binvar x
sol = optimize([-10<=u<=10,implies(x,u>=3+2*x)],-x);
assert(sol.problem == 0)
assert(abs(value(x)-1) <= 1e-4);
assert(value(u)>4.999);

sol = optimize([-10<=u<=10,-6 <= y <= 6,implies(x,[u>=3+2*x, y >= u])],-x);
assert(sol.problem == 0)
assert(abs(value(x)-1) <= 1e-4);
assert(value(u)>4.999);
assert(value(y-u)>-0.0001);

sol = optimize([-10<=u<=10,-4 <= y <= 4,implies(x,[u>=3+2*x, y >= u])],-x);
assert(sol.problem == 0)
assert(abs(value(x)-0) <= 1e-4);
assert(value(u)>=-10);
assert(value(y)>=-4);


binvar x
sdpvar u y
sol = optimize([-10<=u<=10,-4 <= y <= 4,implies(x,[u==3+2*x])],-x);
assert(sol.problem == 0)
assert(abs(value(x)-1) <= 1e-4);
assert(abs(value(u)-5) <= 1e-4);

binvar x
sdpvar u y
sol = optimize([-10<=u<=10,-4 <= y <= 4,implies(x,[u==3+2*x, y >= 2])],-x);
assert(sol.problem == 0)
assert(abs(value(x)-1) <= 1e-4);
assert(abs(value(u)-5) <= 1e-4);
assert(value(y)>=2);

binvar x
sdpvar u y
sol = optimize([-10<=u<=10,-4 <= y <= 4,implies(1-x, u==0.3),implies(x,[u==3+2*x, y >= 6])],-x);
assert(sol.problem == 0)
assert(abs(value(x)-0) <= 1e-4);
assert(abs(value(u)-0.3) <= 1e-4);


sdpvar x u
sol = optimize([-10<=u<=10,-1<=x<=1,implies(x>=0,u==-x),implies(x<=0,u==-2*x)],x);
assert(sol.problem == 0)
assert(abs(value(x)--1) <= 1e-4);
assert(abs(value(u)-2) <= 1e-4);


sol = optimize([-10<=u<=10,-1<=x<=1,implies(x>=0,u==-x),implies(x<=0,u==-2*x)],-x);
assert(sol.problem == 0)
assert(abs(value(x)-1) <= 1e-4);
assert(abs(value(u)--1) <= 1e-4);

sol = optimize([-.5<=u<=.5,-1<=x<=1,implies(x>=0,u==-x),implies(x<=0,u==-2*x)],-x);
assert(sol.problem == 0)
assert(abs(value(x)-.5) <= 1e-4);
assert(abs(value(u)--.5) <= 1e-4);

ops = sdpsettings('quadprog.Algorithm','interior-point-convex');
x = sdpvar(2,1);
sol = optimize([-10<=u<=10,-1<=x<=1,implies(x>=0,u==-x(1)),implies(x<=0,u==-2*x(1))],(x+.1)'*(x+.1),ops)
assert(sol.problem == 0)
assert(abs(value(u)-.2) <= 1e-4);

sol = optimize([-10<=u<=10,-1<=x<=1,implies(x>=0,u==-x(1)),implies(x<=0,u==-2*x(1))],(x-.1)'*(x-.1),ops);
assert(sol.problem == 0)
assert(abs(value(u)--.1) <= 1e-4);

sol = optimize([-10<=u<=10,-1<=x<=1,implies(x>=0,u==-x(1)),implies(x<=0,u==-2*x(1))],(x-[.1;-.1])'*(x-[.1;-.1])+u^2,ops);
assert(sol.problem == 0)
assert(norm(value(x)-[.1;-0.1]) <= 1e-4);
assert(abs(value(u)-0) <= 1e-4);

sdpvar x y u
sol = optimize([-10<=u<=10,-1<=x<=1,-5<=y<=5,implies([x>=0, u==0],y == 4)],(x-1)^2 + u^2,ops);
assert(sol.problem == 0)
assert(abs(value(x)-1) <= 1e-4);
assert(abs(value(u)-0) <= 1e-4);
assert(abs(value(y)-4) <= 1e-4);

binvar d1
sdpvar x y
optimize([-4 <= [x y] <= 10,implies([d1==1,x>=0],y<=1)],-d1-x)
assert(abs(value(d1)-1) <= 1e-4);
assert(value(y) < 1+1e-5)

optimize([-4 <= [x y] <= 10,implies([d1==1,x>=0],y>=1)],-d1-x)
assert(abs(value(d1)-1) <= 1e-4);
assert(value(y) > 1-1e-5)

x = sdpvar(2,1);
u = sdpvar(1);
optimize([-10 <= [x;u] <= 10,implies(x==0.5,u == pi)],(x-.5)'*(x-.5),ops)
assert(norm(value(x)-.5) <= 1e-3);

optimize([-10 <= [x;u] <= 10,implies(x==0.5,u == 15,1e-3)],(x-.5)'*(x-.5),ops)
assert(any(abs(value(x)-0.5)>1e-5));


