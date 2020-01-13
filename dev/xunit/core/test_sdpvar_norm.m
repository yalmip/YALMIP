function tests = test_sdpvar_norm
tests = functiontests(localfunctions);

function test1(dummy)
% To improve performance, we don't introduce auxilliary variables for
% elements which are fixed
sdpvar x y
assign(x,1);
assert(abs(value(norm([x;-3],1)) - 4)<= 1e-8)
optimize(norm([x;-3],1)<=y,y);
assert(norm(value(y) - 3)<= 1e-8)
assert(norm(value(x) - 0)<= 1e-8)
