function tests = test_primalize
tests = functiontests(localfunctions);

function test1(testCase)

n = 10;
randn('seed',123456789);

A = randn(n);A = A - max(real(eig(A)))*eye(n)*1.5;
B = randn(n,1);
C = randn(1,n);

t = sdpvar(1,1);
P = sdpvar(n,n);

obj = t;
F = (kyp(A,B,P,blkdiag(C'*C,-t)) <= 0);

[Fp,objp,free] = primalize(F,-obj);

sol = optimize(Fp,objp,sdpsettings('removeequalities',1,'verbose',0));

testCase.assertTrue(sol.problem == 0)
testCase.assertTrue(abs(value(obj) - 3.38842085719774) <= 1e-5)


function test2(testCase)

n = 10;
randn('seed',123456789);

A = randn(n);A = A - max(real(eig(A)))*eye(n)*1.5;
B = randn(n,1);
C = randn(1,n);

t = sdpvar(1,1);
P = sdpvar(n,n);

obj = t;
F = (kyp(A,B,P,blkdiag(C'*C,-t)) <= 0);

[Fp,objp,free] = primalize(F,-obj);

sol = optimize(Fp,objp,sdpsettings('removeequalities',1,'verbose',0));

testCase.assertTrue(sol.problem == 0)
testCase.assertTrue(abs(value(obj) - 3.38842085719774) <= 1e-5)


function test3(testCase)

n = 50;
randn('seed',123456789);

A = randn(n);A = A - max(real(eig(A)))*eye(n)*1.5;
B = randn(n,1);
C = randn(1,n);

t = sdpvar(1,1);
P = sdpvar(n,n);

obj = t;
F = (kyp(A,B,P,blkdiag(C'*C,-t)) <= 0);

[Fp,objp,free] = primalize(F,-obj);

sol = optimize(Fp,objp,sdpsettings('removeequalities',1,'verbose',0));

testCase.assertTrue(sol.problem == 0)
testCase.assertTrue(abs(value(obj) - 3.74980908287456) <= 1e-5)


function test4(testCase)

n = 50;
randn('seed',123456789);

A = randn(n);A = A - max(real(eig(A)))*eye(n)*1.5;
B = randn(n,1);
C = randn(1,n);

t = sdpvar(1,1);
P = sdpvar(n,n);

obj = t;
F = (kyp(A,B,P,blkdiag(C'*C,-t)) <= 0);

[Fp,objp,free] = primalize(F,-obj);

sol = optimize(Fp,objp,sdpsettings('removeequalities',2,'verbose',0))

testCase.assertTrue(sol.problem == 0)
testCase.assertTrue(abs(value(obj) - 3.74980908287456) <= 1e-5)

