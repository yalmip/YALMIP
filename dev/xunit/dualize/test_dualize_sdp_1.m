function tests = test_sdpvar_dualize_sdp1_1
tests = functiontests(localfunctions);

function test1(dummy)
% TEST 1
A = randn(3,3);A = -A*A';
P = sdpvar(3,3);
t = sdpvar(1,1);
y = sdpvar(1,1);
F = (A'*P+P*A <= -eye(3));
F = F + (P >= A*A') + (P(3,3)>=0) + (t+y >= 7) + (P(2,2)>=4)+(P(1,1:2)>=t) + (t>=12);
obj = trace(P)+y;
    
sol1  = optimize(F,obj);
obj1 = value(obj);
p1   = checkset(F);

sol2 = optimize(F,obj,sdpsettings('dualize',1));
obj2 = value(obj);
p2   = checkset(F);

assert(sol1.problem == 0);
assert(sol2.problem == 0);
assert(abs(obj1 - obj2) <= 1e-4);
assert(abs(min(p1))<= 1e-4)
assert(abs(min(p2))<= 1e-4)