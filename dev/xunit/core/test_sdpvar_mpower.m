function test_sdpvar_power

x = sdpvar(1);
X = sdpvar(2);
Y = randn(2);Y = Y*Y';

assertEqual(x^2-x.^2,0)
assertEqual([x^2 x^2]-[x x].^2,[0 0])
assertEqual(X^2-X*X,zeros(2))
try
    X^2.5;
	assertTrue(false);
catch
end
try
    X^x;
	assertTrue(false);
catch
end
try
    X^X;
	assertTrue(false);
catch
end
try
    X^Y;
	assertTrue(false);
catch
end
try
    Y^X;
	assertTrue(false);
catch
end
assign(x,3);
assertElementsAlmostEqual(double(Y^x),Y^3)
assign(x,1);
assertElementsAlmostEqual(double(Y^-1),Y^-1)

Y = [1 2;3 4];
assign(x,2);
assertElementsAlmostEqual(double(Y^x),Y^2)
