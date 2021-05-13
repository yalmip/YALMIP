function tests = test_ndsdpvar_cat
tests = functiontests(localfunctions);

function test1(testCase)

x = sdpvar(1,2,3);
p = randn(1,2,3);
assign(x,p);

testCase.assertTrue(isequal(value(cat(1,x,x)),cat(1,p,p)))
testCase.assertTrue(isequal(value(cat(2,x,x)),cat(2,p,p)))
testCase.assertTrue(isequal(value(cat(3,x,x)),cat(3,p,p)))

testCase.assertTrue(isequal(value(cat(1,x,x,x)),cat(1,p,p,p)))
testCase.assertTrue(isequal(value(cat(2,x,x,x)),cat(2,p,p,p)))
testCase.assertTrue(isequal(value(cat(3,x,x,x)),cat(3,p,p,p)))

y = sdpvar(1,2,3);
q = randn(1,2,3);
assign(y,q);
testCase.assertTrue(isequal(value(cat(1,x,y)),cat(1,p,q)))
testCase.assertTrue(isequal(value(cat(2,x,y)),cat(2,p,q)))
testCase.assertTrue(isequal(value(cat(3,x,y)),cat(3,p,q)))

testCase.assertTrue(isequal(value(cat(1,p,y)),cat(1,p,q)))
testCase.assertTrue(isequal(value(cat(2,p,y)),cat(2,p,q)))
testCase.assertTrue(isequal(value(cat(3,p,y)),cat(3,p,q)))

testCase.assertTrue(isequal(value(cat(1,x,q)),cat(1,p,q)))
testCase.assertTrue(isequal(value(cat(2,x,q)),cat(2,p,q)))
testCase.assertTrue(isequal(value(cat(3,x,q)),cat(3,p,q)))

z = sdpvar(1,2);
r = randn(1,2);
assign(z,r);
testCase.assertTrue(isequal(value(cat(3,x,z)),cat(3,p,r)))
