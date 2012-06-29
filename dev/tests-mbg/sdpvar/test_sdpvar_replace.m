function test_sdpvar_replace
sdpvar t

p = 1+t^2;
p2 = replace(p,t,2*t);
mbg_asserttrue(isequal(p2-(1+4*t^2),0))
