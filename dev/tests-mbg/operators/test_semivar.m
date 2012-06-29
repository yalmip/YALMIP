function test_abs


A = magic(10);
b = A(:,end);
A = A(:,1:5);
x = semivar(5,1);
e = b-A*(x-1);
obj = norm(e,1);
sol = solvesdp([1 <= x <= 2],obj);

mbg_asserttrue(sol.problem == 0)
mbg_asserttolequal(double(obj), 25.211, 1e-3);

% Still not working
obj = e'*e
sol = solvesdp([1 <= x <= 2],obj);
mbg_asserttrue(sol.problem == 0)
mbg_asserttolequal(double(obj), 133.2742, 1e-3);

