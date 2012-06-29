function test_sdpvar_geomean

x = sdpvar(2,1);
sol = solvesdp(set(0.9 <= x <= 3.3) + set(ismember(x,[1 2 3.2])),x(1)+x(2),sdpsettings('verbose',0));
mbg_asserttolequal(double(x'), [1 1], 1e-5);
sol = solvesdp(set(0.9 <= x <= 3.3) + set(ismember(x,[1 2 3.2])),-x(1)+x(2),sdpsettings('verbose',0));
mbg_asserttolequal(double(x'), [3.2 1], 1e-5);
sol = solvesdp(set(0.9 <= x <= 3.3) + set(ismember(x,[1 2 3.2])),x(1)-x(2),sdpsettings('verbose',0));
mbg_asserttolequal(double(x'), [1 3.2], 1e-5);
sol = solvesdp(set(0.9 <= x <= 3.3) + set(ismember(x,[1 2 3.2])),-x(1)-x(2),sdpsettings('verbose',0));
mbg_asserttolequal(double(x'), [3.2 3.2], 1e-5);







 
