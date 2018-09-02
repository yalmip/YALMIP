function test_cone_in_bmibnb

x = [1 2 3 4 5 6]';
t = (0:0.1:2*pi)';
A = [sin(t) sin(2*t) sin(3*t) sin(4*t) sin(5*t) sin(6*t)];
e = (-4+8*sin(t));
y = A*x+e;
xhat = sdpvar(6,1);
sdpvar u v
F = [-10 <= xhat <= 10, cone(y-A*xhat,u), cone(exp(0.1*xhat.^2),v), exp(xhat) <= exp(4)];
sol = solvesdp(F,u + v,sdpsettings('solver','bmibnb'));

mbg_asserttrue(sol.problem == 0);
mbg_asserttrue(abs(double(u+v)-5.343096e+01) <= 1e-2);
