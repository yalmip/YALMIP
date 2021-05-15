function tests = test_sdpvar_mpower
tests = functiontests(localfunctions);

function test1(testCase)

x = sdpvar(1);
X = sdpvar(2);
Y = randn(2);Y = Y*Y';

testCase.assertTrue(isequal(x^2-x.^2,0))
testCase.assertTrue(isequal([x^2 x^2]-[x x].^2,[0 0]))
testCase.assertTrue(isequal(X^2-X*X,zeros(2)))
try
    X^2.5;
	testCase.assertTrue(false);
catch
end
try
    X^x;
	testCase.assertTrue(false);
catch
end
try
    X^X;
	testCase.assertTrue(false);
catch
end
try
    X^Y;
	testCase.assertTrue(false);
catch
end
try
    Y^X;
	testCase.assertTrue(false);
catch
end
assign(x,3);
testCase.assertTrue(norm(value(Y^x) - Y^3) <= 1e-5)
assign(x,1);
testCase.assertTrue(norm(value(Y^-1) - Y^-1) <= 1e-5)

Y = [1 2;3 4];
assign(x,2);
testCase.assertTrue(norm(value(Y^x) - Y^2) <= 1e-5)
