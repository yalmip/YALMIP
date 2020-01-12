function tests = test_logic_nzz_1
tests = functiontests(localfunctions);

function test1(dummy)

x = sdpvar(4,1);
sol = solvesdp((-10 <= x <= 10) + (nnz(x)==3),(x-2)'*(x-2))
assert(sol.problem == 0);
assert(norm(sort(double(x))-sort([2;2;2;0])) < 1e-3)

sol = solvesdp((-10 <= x <= 10) + (nnz(x)<=3),(x-2)'*(x-2))
assert(sol.problem == 0);
assert(norm(sort(double(x))-sort([2;2;2;0])) < 1e-3)



sdpvar x y z
sol = solvesdp((200>=[x y z] >= -20) + (nnz([x;y;z] == 3) >= 2),2*x+y+z);
assert(sol.problem == 0);
assert(norm(double([x y z] - [-20 3 3])) < 1e-3)
