function tests = test_isconvex
tests = functiontests(localfunctions);

function test1(dummy)

sdpvar x y
assert(isconvex(x+y));
assert(isconvex(x+y^2));
assert(isconvex(exp(x+y)));
assert(isconvex(max(x,exp(x+y))));
assert(~isconvex(-exp(x+y)));
assert(~isconvex(-max(x,exp(x+y))));
assert(isnan(isconvex(max(x,min(x,-x)))));






 
