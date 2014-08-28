function test_sdpvar_geomean_complex

% Test complex matrix geomean, length ~2^n
randn('seed',1234);
rand('seed',1234);
D = randn(5,5)+sqrt(-1)*randn(5,5);D = D + D'+eye(5)*10;
P = sdpvar(5,5,'he','co');
obj = geomean(P);
solvesdp((P <= D),-obj);
assertElementsAlmostEqual(double(obj),9.07516376113709, 'absolute',1e-5);
 
% Test complex matrix geomean, length == 2^n
randn('seed',1234);
rand('seed',1234);
D = randn(8,8)+sqrt(-1)*randn(8,8);D = D + D'+eye(8)*20;
P = sdpvar(8,8,'he','co');
obj = geomean(P);
solvesdp((P <= D),-obj);
assertElementsAlmostEqual(double(obj), 18.42071980565500,'absolute', 1e-5);
 