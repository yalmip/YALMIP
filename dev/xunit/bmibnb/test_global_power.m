function bmibnb_integer

sdpvar x y
ops =sdpsettings('solver','bmibnb');
sol = solvesdp([0 <= [x] <= 5, -5 <= y <= 5, x==.1], (x^y-10)^2,ops)

mbg_asserttrue(sol.problem==0);
mbg_asserttolequal(double(y),-1, 1e-4);