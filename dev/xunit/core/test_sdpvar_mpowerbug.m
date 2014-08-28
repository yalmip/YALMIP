function test_sdpvar_mpowerbug

yalmip('clear')
sdpvar x y
x = x^2;
y = y^2;

p1 = x^2;
p2 = x^2;

% This crashed in R20120806
p3 = (x+y)*(x+y)
assertTrue(true);

