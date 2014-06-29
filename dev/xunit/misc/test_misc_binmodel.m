function test_misc_binmodel

sdpvar y
binvar x

[xy,model] = binmodel(x*y,[2 <= y <= 5]);
solvesdp(model,-xy)
assertElementsAlmostEqual(double(xy), 5, 'absolute', 1e-4);

[xy,model] = binmodel(x*y,abs(y) <= 5);
solvesdp(model,xy)
assertElementsAlmostEqual(double(xy), -5, 'absolute',1e-4);
