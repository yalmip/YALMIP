function tests = test_cone_in_fmincon
tests = functiontests(localfunctions);

function test1(dummy)
x = [1 2 3 4 5 6]';
t = (0:0.1:2*pi)';
A = [sin(t) sin(2*t) sin(3*t) sin(4*t) sin(5*t) sin(6*t)];
e = (-4+8*sin(t));
y = A*x+e;
xhat = sdpvar(6,1);
sdpvar u v
F = [cone(y-A*xhat,u), cone(exp(0.1*xhat.^2),v), exp(xhat) <= exp(4)];
sol = optimize(F,u + v,sdpsettings('solver','fmincon'));

assert(sol.problem == 0);
assert(abs(value(u+v)-5.343096e+01) <= 1e-3);
