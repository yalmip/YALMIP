function tests = test_sdpvar_norm
tests = functiontests(localfunctions);

function test1(dummy)
% To improve performance, we don't introduce auxilliary variables for
% elements which are fixed
sdpvar x y
assign(x,1);
assert(abs(double(norm([x;-3],1)) - 4)<= 1e-8)
solvesdp(norm([x;-3],1)<=y,y);
assert(norm(double(y) - 3)<= 1e-8)
assert(norm(double(x) - 0)<= 1e-8)
