function test_isconvex

sdpvar x y
mbg_asserttrue(isconvex(x+y));
mbg_asserttrue(isconvex(x+y^2));
mbg_asserttrue(isconvex(exp(x+y)));
mbg_asserttrue(isconvex(max(x,exp(x+y))));
mbg_asserttrue(~isconvex(-exp(x+y)));
mbg_asserttrue(~isconvex(-max(x,exp(x+y))));
mbg_asserttrue(isnan(isconvex(max(x,min(x,-x)))));






 
