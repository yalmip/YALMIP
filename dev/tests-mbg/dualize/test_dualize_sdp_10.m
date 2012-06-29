function sdp5

K = sdpvar(4);
sol = solvesdp([K>=0,K(2,3)==1],trace(K)+norm(K(:)),sdpsettings('dualize',1));

mbg_asserttolequal(sol.problem,0);
