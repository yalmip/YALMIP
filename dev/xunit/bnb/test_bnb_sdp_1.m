function test_bnb_sdp_1

% Test two regression bugs
% 1. quadratic costs in binary SDP
% 2. OR on SDP constraints

X = sdpvar(3,3);
obj = trace((X-2*eye(3))*(X-2*eye(3))');
sol = solvesdp((X<=3*eye(3)) + ((X>=eye(3)) | (X<=-eye(3))) + (-50 <= X(:) <= 50),obj)

mbg_asserttrue(sol.problem==0);
mbg_asserttolequal(double(obj), 0, 1e-5);

