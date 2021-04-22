function tests = test_mpt_norm_repeated
tests = functiontests(localfunctions);

function test1(dummy)

% As simple as it gets...
yalmip('clear')
sdpvar u0
x0 = sdpvar(2,1);

bounds(x0,-2,2);
bounds(u0,-2,2);
x1 = [0.1 0.2;0.3 0.4]*x0 + [0.5;0.6]*u0
F = ([]);
F = F + (norm(x1,1) <= norm(x0,1));
F = F +  (-2<=x0<=2);
F = F +  (-2<=x1<=2);
F = F +  (-2<=u0<=2);
obj = norm(x0,1);
[sol,dgn,VV,JJ] = solvemp(F, obj, [], x0)

assert(length(sol) == 4)
assert(dgn.problem == 0)

assign(x0,[1;1])
assert(abs(value(JJ)-2) <= 1e-5)

yalmip('clear')
sdpvar u0 u1
x0 = sdpvar(2,1);
x1 = sdpvar(2,1);
x2 = sdpvar(2,1);
bounds(x0,-2,2);
bounds(x1,-2,2);
bounds(x2,-2,2);
bounds(u0,-2,2);
bounds(u1,-2,2);
%x1 = [0.1 0.2;0.3 0.4]*x0 + [0.5;0.6]*u0
F = (x1 ==[0.1 0.2;0.3 0.4]*x0 + [0.5;0.6]*u0);
F = F + (x2 ==[0.1 0.2;0.3 0.4]*x1 + [0.5;0.6]*u1);
F = F + (norm(x1,1) <= norm(x0,1));
F = F + (norm(x2,1) <= norm(x1,1));
F = F +  (-2<=x0<=2);
F = F +  (-2<=x1<=2);
F = F +  (-2<=x2<=2);
F = F +  (-2<=u0<=2);
obj = norm(x0,1)+norm(x1,1);
sol = solvemp(F, obj, [], x0)
assert(length(sol) == 6);

yalmip('clear')
sdpvar u0 u1
x0 = sdpvar(2,1);
x1 = sdpvar(2,1);
x2 = sdpvar(2,1);
bounds(x0,-2,2);
bounds(x1,-2,2);
bounds(x2,-2,2);
bounds(u0,-2,2);
bounds(u1,-2,2);
F = (x1 ==[0.1 0.2;0.3 0.4]*x0 + [0.5;0.6]*u0);
F = F + (x2 ==[0.1 0.2;0.3 0.4]*x1 + [0.5;0.6]*u1);
F = F + (norm(x2,1) <= norm(x1,1));
F = F + (norm(x1,1) <= norm(x0,1));
F = F +  (-2<=x0<=2);
F = F +  (-2<=x1<=2);
F = F +  (-2<=x2<=2);
F = F +  (-2<=u0<=2);
obj = norm(x0,1)+norm(x1,1);
sol = solvemp(F, obj, [], x0)
assert(length(sol) == 6)