function tests = test_sdpvar_mpower
tests = functiontests(localfunctions);

function test1(dummy)

x = sdpvar(1);
X = sdpvar(2);
Y = randn(2);Y = Y*Y';

assert(isequal(x^2-x.^2,0))
assert(isequal([x^2 x^2]-[x x].^2,[0 0]))
assert(isequal(X^2-X*X,zeros(2)))
try
    X^2.5;
	assert(false);
catch
end
try
    X^x;
	assert(false);
catch
end
try
    X^X;
	assert(false);
catch
end
try
    X^Y;
	assert(false);
catch
end
try
    Y^X;
	assert(false);
catch
end
assign(x,3);
assert(norm(value(Y^x) - Y^3) <= 1e-5)
assign(x,1);
assert(norm(value(Y^-1) - Y^-1) <= 1e-5)

Y = [1 2;3 4];
assign(x,2);
assert(norm(value(Y^x) - Y^2) <= 1e-5)
