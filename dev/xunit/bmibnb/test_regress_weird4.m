function test_regress_weird4

% Issue #188
yalmip('clear')
sdpvar x y 
p = sin(1+y*x)^2+cos(y*x);
sdpvar z w
sol = solvesdp([-1 <= [x y z w] <= 1, p <= 3],p,sdpsettings('solver','bmibnb'));

mbg_asserttolequal(sol.problem,0);
mbg_asserttolequal(double(p),.5403, 1e-2);