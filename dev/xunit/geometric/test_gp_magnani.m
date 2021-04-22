function tests = test_gp_magnani
tests = functiontests(localfunctions);

function test1(dummy)

x=sdpvar(1,1);
y=sdpvar(1,1);
t=x/y;
F=(t>=0.5);
F=F+(y<=2);
F=F+(y>=1);
obj=x^2*y^3;
F=F+(t<=1);
F=F+(t>=1);
F=F+(y^2<=4);
F = [F, x>=0, y>=0];
sol = optimize(F,obj);
assert(sol.problem == 0);
assert(abs(value(obj)-1)<=1e-4);
assert(all(abs(value([x y t]) - [1 1 1]) <= 1e-4));

sol = optimize(F,1/obj);
assert(sol.problem == 0);
assert(abs(value(obj)-32) <= 1e-3);
assert(all(abs(value([x y t]) - [2 2 1]) <= 1e-4));

function test2(dummy)

x=sdpvar(1,1);
y=sdpvar(1,1);
t=x/y;
F=(t>=0.5);
F=F+(y<=2);
F=F+(y>=1);
obj=x^2*y^3;
F=F+(t<=1);
F=F+(t>=1);
F=F+(y^2.5<=4);
F = [F, x>=0];
sol = optimize(F,obj);
assert(sol.problem == 0);
assert(abs(value(obj)-1)<=1e-4);
assert(all(abs(value([x y t]) - [1 1 1]) <= 1e-4));

sol = optimize(F,1/obj);
assert(sol.problem == 0);
assert(abs(value(obj)-16) <= 1e-3);
assert(all(abs(value([x y t]) - [1.74110112659225   1.74110112659225   1.00000000000000]) <=  1e-4));