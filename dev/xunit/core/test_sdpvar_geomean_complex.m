function tests = test_geomean_complex
tests = functiontests(localfunctions);

function test1(dummy)
% Test complex matrix geomean, length ~2^n
randn('seed',1234);
rand('seed',1234);
D = randn(5,5)+sqrt(-1)*randn(5,5);D = D + D'+eye(5)*10;
P = sdpvar(5,5,'he','co');
obj = geomean(P);
optimize((P <= D),-obj);
assert(abs(value(obj) - 9.07516376113709) <= 1e-5)
 
% Test complex matrix geomean, length == 2^n
randn('seed',1234);
rand('seed',1234);
D = randn(8,8)+sqrt(-1)*randn(8,8);D = D + D'+eye(8)*20;
P = sdpvar(8,8,'he','co');
obj = geomean(P);
optimize((P <= D),-obj);
assert((value(obj) - 18.42071980565500) <= 1e-5)
 