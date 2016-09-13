function test_operator_interp1

yalmip('clear');
sdpvar x
xv = -1:.1:1;
y = xv.^2+xv+sin(xv*3)*3;
optimize([],interp1(xv,y,x),sdpsettings('solver','bmibnb'))