function tests = test_operator_interp1
tests = functiontests(localfunctions);

function test1(dummy)
yalmip('clear');
sdpvar x
xv = -1:.2:1;
y = xv.^2+xv+sin(xv*3)*3;
sol = optimize([],interp1(xv,y,x),sdpsettings('solver','bmibnb'));
assert(sol.problem == 0);
assert(abs(value(x)--.6)<1e-3);
sol = optimize([],interp1(xv,y,x,'spline'),sdpsettings('solver','bmibnb','bmibnb.relgaptol',1e-6,'bmibnb.absgaptol',1e-4));
assert(sol.problem == 0);
assert(abs(value(x)--.52059)<1e-2);