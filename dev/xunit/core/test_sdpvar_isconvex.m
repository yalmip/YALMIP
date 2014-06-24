function test_isconvex

sdpvar x y
assertTrue(isconvex(x+y));
assertTrue(isconvex(x+y^2));
assertTrue(isconvex(exp(x+y)));
assertTrue(isconvex(max(x,exp(x+y))));
assertTrue(~isconvex(-exp(x+y)));
assertTrue(~isconvex(-max(x,exp(x+y))));
assertTrue(isnan(isconvex(max(x,min(x,-x)))));






 
