function test_operator_entropy

x = sdpvar(2,1);
assign(x,1);
sol = solvesdp((x >= 0.1),-entropy(x),sdpsettings('usex0',1));

assertTrue(sol.problem == 0);
assertElementsAlmostEqual(double(sum(x)),0.73,'absolute',1e-2);
assertElementsAlmostEqual(-entropy(double(x)),-0.7357588,'absolute',1e-3);

x1 = sdpvar(1,1);
y = sdpvar(1,1);
x2 = sdpvar(1,1);
x = [x1;x2];

assign(x,1);
sol = solvesdp(([x;y] >= 0.1),-entropy(x)+y,sdpsettings('usex0',1));

assertTrue(sol.problem == 0);
assertElementsAlmostEqual(double(sum(x)),0.73,'absolute',1e-2);
