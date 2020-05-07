function tests = test_bnb_sdp_1
tests = functiontests(localfunctions);

function test1(dummy)
% Test two regression bugs
% 1. quadratic costs in binary SDP
% 2. OR on SDP constraints

X = sdpvar(3,3);
obj = trace((X-2*eye(3))*(X-2*eye(3))');
sol = optimize((X<=3*eye(3)) + ((X>=eye(3)) | (X<=-eye(3))) + (-50 <= X(:) <= 50),obj)

assert(sol.problem==0);
assert(abs(value(obj)) <= 1e-5);

