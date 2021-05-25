function tests = test_slayer_linear
tests = functiontests(localfunctions);

function test_lp_indisguise(testCase)

sdpvar x y
options = sdpsettings('solver','fmincon','debug',1);

P = (ones(4)+eye(4)*3)/3;
F = P*diag([x+1;1-y;1-x;y+1])*P>=0;
options.slayer.algorithm = 0;
sol = optimize(F,-x-y,options);
testCase.assertTrue(sol.problem == 0)
options.slayer.algorithm = 1;
sol = optimize(F,-x-y,options)
testCase.assertTrue(sol.problem == 0)

F = diag([x+1;1-y;1-x;y+1])>=0;
options.slayer.algorithm = 0;
sol = optimize(F,-x-y,options);
testCase.assertTrue(sol.problem == 0)
options.slayer.algorithm = 1;
sol = optimize(F,-x-y,options)
testCase.assertTrue(sol.problem == 0)

A = [1 0;0.4 1];
B = [0.4;0.08];
L = [1.9034 1.1501];
Y = sdpvar(2,2);
F = [Y Y*(A-B*L)';(A-B*L)*Y Y]>=0;
F = F+[L*Y*L'<=1];
options.slayer.algorithm = 0;
sol = optimize(F,-trace(Y),options);
testCase.assertTrue(sol.problem == 0)
testCase.assertTrue(abs(trace(value(Y))-10.26038)<=1e-3)
options.slayer.algorithm = 1;
sol = optimize(F,-trace(Y),options);
testCase.assertTrue(sol.problem == 0)
testCase.assertTrue(abs(trace(value(Y))-10.26038)<=1e-3)

W=[Y Y*(A-B*L)';(A-B*L)*Y Y];
F = blkdiag(W,W,Y)>=0;
F = F+[L*Y*L'<=1] + [blkdiag(Y,2*Y)>=0]+[10*Y>=0];
options.slayer.algorithm = 0;
sol = optimize(F,-trace(Y),options);
testCase.assertTrue(sol.problem == 0)
testCase.assertTrue(abs(trace(value(Y))-10.26038)<=1e-3)
options.slayer.algorithm = 1;
sol = optimize(F,-trace(Y),options);
testCase.assertTrue(sol.problem == 0)
testCase.assertTrue(abs(trace(value(Y))-10.26038)<=1e-3)

options.slayer.algorithm = 0;
sol = optimize(F,-logdet(Y),options);
testCase.assertTrue(sol.problem == 0)
testCase.assertTrue(abs(value(geomean(Y))-1.27059)<=1e-3)
options.slayer.algorithm = 1;
sol = optimize(F,-logdet(Y),options);
testCase.assertTrue(sol.problem == 0)
testCase.assertTrue(abs(value(geomean(Y))-1.27059)<=1e-3)

