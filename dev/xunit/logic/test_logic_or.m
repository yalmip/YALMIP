function tests = test_logic_or
tests = functiontests(localfunctions);

function test1(dummy)
binvar y1 y2 y3
optimize([true(or(y1,y2,y3))],y1+y2+y3)
assert(sum(value([y1 y2 y3]))==1);
