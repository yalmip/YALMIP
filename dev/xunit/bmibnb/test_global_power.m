function test_global_power

sdpvar x y
ops =sdpsettings('solver','bmibnb');
sol = solvesdp([0 <= [x] <= 5, -5 <= y <= 5, x==.1], (x^y-10)^2,ops)

mbg_asserttrue(sol.problem==0);
mbg_asserttolequal(double(y),-1, 1e-4);

x = sdpvar(2,1);
y = sdpvar(2,1);
ops =sdpsettings('solver','bmibnb');
sol = solvesdp([0 <= [x] <= 5, -5 <= y <= 5, x==.1], sum(sum((x.^y-10).^2)),ops)

mbg_asserttrue(sol.problem==0);
mbg_asserttolequal(double(y),[-1;-1], 1e-4);