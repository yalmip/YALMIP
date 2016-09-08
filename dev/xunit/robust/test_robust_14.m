function test_robust_14

sdpvar t x 
w = sdpvar(2,1);
M = momentmodel([w'*w == 1, w>=0.5])
P = sosmodel(sos(1 + w(1)*x + w(2)*x^2 - t),[],[],[t;w]);
sol = optimize([M,P,uncertain(M)],-t)
mbg_asserttrue(sol.problem == 0);
mbg_asserttolequal(double(t),.625, 1e-3);

