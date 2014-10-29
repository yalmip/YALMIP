function test_operator_implies

% Binary variable implies LP constriants
sdpvar y u
binvar x
sol = solvesdp([-10<=u<=10,implies(x,u>=3+2*x)],-x);
assertTrue(sol.problem == 0)
mbg_asserttolequal(double(x),1,1e-4);
assertTrue(double(u)>4.999);

sol = solvesdp([-10<=u<=10,-6 <= y <= 6,implies(x,[u>=3+2*x, y >= u])],-x);
assertTrue(sol.problem == 0)
mbg_asserttolequal(double(x),1,1e-4);
assertTrue(double(u)>4.999);
assertTrue(double(y-u)>-0.0001);

sol = solvesdp([-10<=u<=10,-4 <= y <= 4,implies(x,[u>=3+2*x, y >= u])],-x);
assertTrue(sol.problem == 0)
mbg_asserttolequal(double(x),0,1e-4);
assertTrue(double(u)>=-10);
assertTrue(double(y)>=-4);


binvar x
sdpvar u y
sol = solvesdp([-10<=u<=10,-4 <= y <= 4,implies(x,[u==3+2*x])],-x);
assertTrue(sol.problem == 0)
mbg_asserttolequal(double(x),1,1e-4);
mbg_asserttolequal(double(u),5,1e-4);

binvar x
sdpvar u y
sol = solvesdp([-10<=u<=10,-4 <= y <= 4,implies(x,[u==3+2*x, y >= 2])],-x);
assertTrue(sol.problem == 0)
mbg_asserttolequal(double(x),1,1e-4);
mbg_asserttolequal(double(u),5,1e-4);
mbg_asserttrue(double(y)>=2);

binvar x
sdpvar u y
sol = solvesdp([-10<=u<=10,-4 <= y <= 4,implies(1-x, u==0.3),implies(x,[u==3+2*x, y >= 6])],-x);
assertTrue(sol.problem == 0)
mbg_asserttolequal(double(x),0,1e-4);
mbg_asserttolequal(double(u),0.3,1e-4);


sdpvar x u
sol = solvesdp([-10<=u<=10,-1<=x<=1,implies(x>=0,u==-x),implies(x<=0,u==-2*x)],x);
assertTrue(sol.problem == 0)
mbg_asserttolequal(double(x),-1,1e-4);
mbg_asserttolequal(double(u), 2,1e-4);


sol = solvesdp([-10<=u<=10,-1<=x<=1,implies(x>=0,u==-x),implies(x<=0,u==-2*x)],-x);
assertTrue(sol.problem == 0)
mbg_asserttolequal(double(x),1,1e-4);
mbg_asserttolequal(double(u),-1,1e-4);

sol = solvesdp([-.5<=u<=.5,-1<=x<=1,implies(x>=0,u==-x),implies(x<=0,u==-2*x)],-x);
assertTrue(sol.problem == 0)
mbg_asserttolequal(double(x),0.5,1e-4);
mbg_asserttolequal(double(u),-0.5,1e-4);

ops = sdpsettings('quadprog.Algorithm','interior-point-convex');
x = sdpvar(2,1);
sol = solvesdp([-10<=u<=10,-1<=x<=1,implies(x>=0,u==-x(1)),implies(x<=0,u==-2*x(1))],(x+.1)'*(x+.1),ops)
assertTrue(sol.problem == 0)
mbg_asserttolequal(double(u),0.2,1e-3);

sol = solvesdp([-10<=u<=10,-1<=x<=1,implies(x>=0,u==-x(1)),implies(x<=0,u==-2*x(1))],(x-.1)'*(x-.1),ops);
assertTrue(sol.problem == 0)
mbg_asserttolequal(double(u),-.1,1e-3);


sol = solvesdp([-10<=u<=10,-1<=x<=1,implies(x>=0,u==-x(1)),implies(x<=0,u==-2*x(1))],(x-[.1;-.1])'*(x-[.1;-.1])+u^2,ops);
assertTrue(sol.problem == 0)
mbg_asserttolequal(double(x),[.1;-0.1],1e-3);
mbg_asserttolequal(double(u),0,1e-3);

sdpvar x y u
sol = solvesdp([-10<=u<=10,-1<=x<=1,-5<=y<=5,implies([x>=0, u==0],y == 4)],(x-1)^2 + u^2,ops);
assertTrue(sol.problem == 0)
mbg_asserttolequal(double(x),1,2e-3);
mbg_asserttolequal(double(u),0,2e-3);
mbg_asserttolequal(double(y),4,2e-3);

binvar d1
sdpvar x y
solvesdp([-4 <= [x y] <= 10,implies([d1==1,x>=0],y<=1)],-d1-x)
mbg_asserttolequal(double(d1),1,1e-4);
assertTrue(double(y) < 1+1e-5)

solvesdp([-4 <= [x y] <= 10,implies([d1==1,x>=0],y>=1)],-d1-x)
mbg_asserttolequal(double(d1),1,1e-4);
assertTrue(double(y) > 1-1e-5)

x = sdpvar(2,1);
u = sdpvar(1);
solvesdp([-10 <= [x;u] <= 10,implies(x==0.5,u == pi)],(x-.5)'*(x-.5),ops)
mbg_asserttolequal(norm(double(x)-.5),0,1e-3);

solvesdp([-10 <= [x;u] <= 10,implies(x==0.5,u == 15,1e-3)],(x-.5)'*(x-.5),ops)
mbg_asserttrue(any(abs(double(x)-0.5)>1e-5));


