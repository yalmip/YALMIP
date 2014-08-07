function test_regress_weird5

% Issue #187
yalmip('clear')
sdpvar x y
p = sin(1+y)^2+cos(y*x);
sol = solvesdp([-1 <= [x y] <= 1, p <= 3],p,sdpsettings('solver','bmibnb'));
mbg_asserttolequal(sol.problem,0);
mbg_asserttolequal(double(p),.5403, 1e-2);