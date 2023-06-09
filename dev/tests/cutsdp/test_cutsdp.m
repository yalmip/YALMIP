function tests = test_cutsdp
tests = functiontests(localfunctions);

function test_diw(testCase)
[F,h] = loadsdpafile('diw_15.dat-s');
sol = optimize(F,h,sdpsettings('solver','cutsdp'));
testCase.assertTrue(sol.problem == 0);
testCase.assertTrue(abs(value(h) - -95) <= 1e-2);
function test_colon(testCase)
[F,h] = loadsdpafile('coloncancer_1_100_5.dat-s')
sol = optimize(F,h,sdpsettings('solver','cutsdp'));
testCase.assertTrue(sol.problem == 0);
testCase.assertTrue(abs(value(h) - 127.47) <= 1e-2);
function test_clique(testCase)
[F,h] = loadsdpafile('clique_20_k3_6_7.dat-s')
sol = optimize(F,h,sdpsettings('solver','cutsdp'));
testCase.assertTrue(sol.problem == 0);
testCase.assertTrue(abs(value(h) - 147) <= 1e-2);
function test_bar(testCase)
[F,h] = loadsdpafile('2x3_3bars.dat-s')
sol = optimize(F,h,sdpsettings('solver','cutsdp'));
testCase.assertTrue(sol.problem == 0);
testCase.assertTrue(abs(value(h) - 2.118) <= 1e-2);