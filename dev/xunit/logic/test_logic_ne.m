function tests = test_logic_ne
tests = functiontests(localfunctions);

function test1(dummy)
yalmip('clear')

binvar x y
optimize([0 <= [x,y] <= 1, x+y~=1],(x+y-1)^2-10*x)
assert(abs(value(x)-1) <= 1e-4)
assert(abs(value(y)-1) <= 1e-4)

binvar x y
optimize([0 <= [x,y] <= 1, x+y~=1],(x+y-1)^2)
assert(abs(value((x+y-1)^2)-1) <= 1e-4)

binvar x1 y1 x2 y2
optimize([0 <= [x1,y1,x2,y2] <= 1, [x1;x2]'+[y1;y2]'~=1],(x1+y1-1)^2-10*x1+(x2+y2-1)^2-10*x2)
assert(abs(value(x1)-1) <= 1e-4)
assert(abs(value(y1)-1) <= 1e-4)