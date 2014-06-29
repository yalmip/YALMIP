function test_ndsdpvar_cat

x = sdpvar(1,2,3);
p = randn(1,2,3);
assign(x,p)
