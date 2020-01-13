function tests = test_sdpvar_dualize_sdp_10
tests = functiontests(localfunctions);

function test1(dummy)

K = sdpvar(4);
sol = optimize([K>=0,K(2,3)==1],trace(K)+norm(K(:)),sdpsettings('dualize',1));

assert(sol.problem == 0);
