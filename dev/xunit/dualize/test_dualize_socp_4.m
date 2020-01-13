function tests = test_sdpvar_dualize_socp_4
tests = functiontests(localfunctions);

function test1(dummy)

X = sdpvar(3,3);
x = sdpvar(3,1);
obj = trace(X)+sum(x);
F = (X>=0) + (cone(1-x(2:end),1+x(1))) + (trace(X)==x(1)+2*x(2)+3*x(3)+4)+(X(1,3)==8);

sol1  = optimize(F,obj);
obj1 = value(obj);
p1   = checkset(F);

sol2 = optimize(F,obj,sdpsettings('dualize',1));
obj2 = value(obj);
p2   = checkset(F);

assert(sol1.problem == 0);
assert(sol2.problem == 0);
assert(abs(obj1 - obj2) <= 1e-4);
assert(abs(min(p1))<= 1e-4)
assert(abs(min(p2))<= 1e-4)