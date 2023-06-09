function tests = test_diff
tests = functiontests(localfunctions);

function test1(testCase)

testCase.assertTrue(isempty(diff(sdpvar(1,1))));
U = sdpvar(1,1);
testCase.assertTrue(isequal(0,diff([U U])));
testCase.assertTrue(isequal(0,diff([U;U])));
testCase.assertTrue(isequal(U - diff([0 U]),0));
testCase.assertTrue(isequal(-U - diff([U 0]),0));
U = sdpvar(2,4);
r = randn(2,4);
assign(U,r)
testCase.assertTrue(norm(value(diff(U)-diff(r)))<1e-12)
testCase.assertTrue(norm(value(diff(U,1)-diff(r,1)))<1e-12)
testCase.assertTrue(norm(value(diff(U,2)-diff(r,2)))<1e-12)
testCase.assertTrue(norm(value(diff(U,1,2)-diff(r,1,2)))<1e-12)
testCase.assertTrue(norm(value(diff(U,2,2)-diff(r,2,2)))<1e-12)
U = sdpvar(4,2);
r = randn(4,2);
assign(U,r)
testCase.assertTrue(norm(value(diff(U)-diff(r)))<1e-12)
testCase.assertTrue(norm(value(diff(U,1)-diff(r,1)))<1e-12)
testCase.assertTrue(norm(value(diff(U,2)-diff(r,2)))<1e-12)
testCase.assertTrue(norm(value(diff(U,1,2)-diff(r,1,2)))<1e-12)
testCase.assertTrue(norm(value(diff(U,2,2)-diff(r,2,2)))<1e-12)
U = sdpvar(4,4);
r = randn(4,4);r=r+r';
assign(U,r)
testCase.assertTrue(norm(value(diff(U)-diff(r)))<1e-12)
testCase.assertTrue(norm(value(diff(U,1)-diff(r,1)))<1e-12)
testCase.assertTrue(norm(value(diff(U,2)-diff(r,2)))<1e-12)
testCase.assertTrue(norm(value(diff(U,1,2)-diff(r,1,2)))<1e-12)
testCase.assertTrue(norm(value(diff(U,2,2)-diff(r,2,2)))<1e-12)
U = sdpvar(4,1);
r = randn(4,1);
assign(U,r)
testCase.assertTrue(norm(value(diff(U)-diff(r)))<1e-12)
testCase.assertTrue(norm(value(diff(U,1)-diff(r,1)))<1e-12)
testCase.assertTrue(norm(value(diff(U,2)-diff(r,2)))<1e-12)
testCase.assertTrue(norm(value(diff(U,1,2)-diff(r,1,2)))<1e-12)
testCase.assertTrue(norm(value(diff(U,2,2)-diff(r,2,2)))<1e-12)
U = sdpvar(4,3,'fu','co');
r = randn(4,3)+randn(4,3)*sqrt(-1);
assign(U,r)
testCase.assertTrue(norm(value(diff(U)-diff(r)))<1e-12)
testCase.assertTrue(norm(value(diff(U,1)-diff(r,1)))<1e-12)
testCase.assertTrue(norm(value(diff(U,2)-diff(r,2)))<1e-12)
testCase.assertTrue(norm(value(diff(U,1,2)-diff(r,1,2)))<1e-12)
testCase.assertTrue(norm(value(diff(U,2,2)-diff(r,2,2)))<1e-12)