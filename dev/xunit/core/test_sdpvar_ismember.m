function test_sdpvar_ismember

x = sdpvar(2,1);
sol = solvesdp((0.9 <= x <= 3.3) + (ismember(x,[1 2 3.2])),x(1)+x(2),sdpsettings('verbose',0));
assertElementsAlmostEqual(double(x'), [1 1], 'absolute',1e-5);
sol = solvesdp((0.9 <= x <= 3.3) + (ismember(x,[1 2 3.2])),-x(1)+x(2),sdpsettings('verbose',0));
assertElementsAlmostEqual(double(x'), [3.2 1],'absolute', 1e-5);
sol = solvesdp((0.9 <= x <= 3.3) + (ismember(x,[1 2 3.2])),x(1)-x(2),sdpsettings('verbose',0));
assertElementsAlmostEqual(double(x'), [1 3.2],'absolute', 1e-5);
sol = solvesdp((0.9 <= x <= 3.3) + (ismember(x,[1 2 3.2])),-x(1)-x(2),sdpsettings('verbose',0));
assertElementsAlmostEqual(double(x'), [3.2 3.2],'absolute', 1e-5);







 
