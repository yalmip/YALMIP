function test_sdpvar_replace

sdpvar t
p = 1+t^2;
p2 = replace(p,t,2*t);
assertTrue(isequal(p2-(1+4*t^2),0))

% Checks that the 0^0 bug in MATLAB6.5 LINUX
% is avoided
sdpvar x t
p = x^2+t;
y = replace(p,t,0);
assertEqual(getbase(y), getbase(x^2));