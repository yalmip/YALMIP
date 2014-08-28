function test_operator_sort

x = sdpvar(4,1);
z = sdpvar(4,1);

[y,loc] = sort(x);

w = randn(4,1);

sol = solvesdp((-100 <= x <= 100)+(z == y),norm(x-w,1));

assertTrue(sol.problem == 0);
assertElementsAlmostEqual(norm(sort(w)-double(z)),0,'absolute',1e-4);


A = ones(20,5);
b = (1:20)';
x = sdpvar(5,1);
e = b-A*x;
F = (mean(x) == median(x)) + (-100 <= x <= 100);
sol = solvesdp(F,norm(e,1));
assertTrue(sol.problem == 0);
assertTrue(abs(mean(double(x))-median(double(x)))<1e-4);
assertTrue(abs(norm(double(e),1)-100)<1e-4);
