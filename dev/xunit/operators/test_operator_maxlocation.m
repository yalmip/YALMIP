function test_operator_maxlocation

% Test for feature #313
y = (1:10)';
x = intvar(10,1);
[val,loc] = max(x);
sol = optimize([sort(x) == y, loc == 5],x'*x)
assertTrue(sol.problem == 0)
assertElementsAlmostEqual(value(x(5)),10,'absolute', 1e-4);

[val,loc] = min(x);
sol = optimize([sort(x) == y, loc == 5],x'*x)
assertTrue(sol.problem == 0)
assertElementsAlmostEqual(value(x(5)),1,'absolute', 1e-4);
