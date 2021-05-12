function tests = test_sdpvar_power
tests = functiontests(localfunctions);

function test1(dummy)
x = sdpvar(1);
y = ([1 1 1]*x).^(0:2);
assert(isequal(getbase(y), getbase([1 x x^2])))
x = sdpvar(1);
y = x.^(0:2);
assert(isequal(getbase(y), getbase([1 x x^2])))

assert(isequal([x 1].^[0 x]-[1 1],zeros(1,2)))
assert(isequal([x 1].^[1 x]-[x 1],zeros(1,2)))
assert(isequal([1].^x-[1],0))
assert(isequal(2.^x-2.^x,0))
assert(isequal(2.^[x x]-[2 2].^x,[0 0]))
assert(isequal([2 3].^[x x]-[2^x 3^x],[0 0]))
assert(isequal([2 3;4 5].^[x x;x x]-[2^x 3^x;4^x 5^x],[0 0;0 0]))

try
    [1 x].^[1 x 1]
	assert(false);
catch
end