function test_sdpvar_times

sdpvar x y
assertEqual(double(isequal([x*x x*x;x*x x*x]-x.*[x x;x x],zeros(2))),1)
assertEqual(double(isequal([x*x x*x;x*x x*x]-[x x;x x].*x,zeros(2))),1)
assertEqual(double(isequal([x*y x*y;x*y x*y]-[x x;x x].*y,zeros(2))),1)
assertEqual(double(isequal([x*y x*y;x*y x*y]-y.*[x x;x x],zeros(2))),1)

x = sdpvar(2,1,'full','complex');
w = randn(2,1) + sqrt(-1)*randn(2,1);
v = randn(2,1) + sqrt(-1)*randn(2,1);
assign(x,w);
assertTrue(all(double(w.*v - x.*v) < 1e-12));