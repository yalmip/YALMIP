function test_sort

x = sdpvar(2,1);

assign(x,1);
sol = solvesdp(set(x >= 0.1),-entropy(x),sdpsettings('usex0',1));

mbg_asserttrue(sol.problem == 0);
mbg_asserttolequal(double(sum(x)),0.73,1e-2);
mbg_asserttolequal(-entropy(double(x)),-0.7357588,1e-3);

x1 = sdpvar(1,1);
y = sdpvar(1,1);
x2 = sdpvar(1,1);
x = [x1;x2];

assign(x,1);
sol = solvesdp(set([x;y] >= 0.1),-entropy(x)+y,sdpsettings('usex0',1));

mbg_asserttrue(sol.problem == 0);
mbg_asserttolequal(double(sum(x)),0.73,1e-2);
