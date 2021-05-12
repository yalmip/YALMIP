function tests = test_sdpvar_ismember
tests = functiontests(localfunctions);

function test1(dummy)
x = sdpvar(2,1);
sol = optimize((0.9 <= x <= 3.3) + (ismember(x,[1 2 3.2])),x(1)+x(2),sdpsettings('verbose',0));
assert(norm(value(x') - [1 1]) <= 1e-5)
sol = optimize((0.9 <= x <= 3.3) + (ismember(x,[1 2 3.2])),-x(1)+x(2),sdpsettings('verbose',0));
assert(norm(value(x') - [3.2 1]) <= 1e-5)
sol = optimize((0.9 <= x <= 3.3) + (ismember(x,[1 2 3.2])),x(1)-x(2),sdpsettings('verbose',0));
assert(norm(value(x') - [1 3.2]) <= 1e-5)
sol = optimize((0.9 <= x <= 3.3) + (ismember(x,[1 2 3.2])),-x(1)-x(2),sdpsettings('verbose',0));
assert(norm(value(x') - [3.2 3.2]) <= 1e-5)







 
