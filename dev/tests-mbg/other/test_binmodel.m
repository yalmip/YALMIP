function test_binmodel

sdpvar y
binvar x

[xy,model] = binmodel(x*y,[2 <= y <= 5]);
solvesdp(model,-xy)
mbg_asserttolequal(double(xy), 5, 1e-4);

[xy,model] = binmodel(x*y,abs(y) <= 5);
solvesdp(model,xy)
mbg_asserttolequal(double(xy), -5, 1e-4);
