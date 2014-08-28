function test_dualize_sdp_8

sdpvar t
A = randn(3,3);A = -A*A';
P = sdpvar(3,3);
%t = sdpvar(1,1);
y = sdpvar(1,1);
F = (A'*P+P*A <= -eye(3));
F = F + (P >= A*A') + (P(3,3)>=0) + (P>=0) + (t+y >= 7) + (t+y >= 7) + (P(2,2)>=4)+(P(1,1:2)>=t) + (t>=12)+(t>=-12);
obj = trace(P)+y+t;

sol1  = solvesdp(F,obj);
obj1 = double(obj);
p1   = checkset(F);

sol2 = solvesdp(F,obj,sdpsettings('dualize',1));
obj2 = double(obj);
p2   = checkset(F);

mbg_asserttolequal(sol1.problem,0);
mbg_asserttolequal(sol2.problem,0);
mbg_asserttolequal(obj1,obj2, 1e-4);
mbg_asserttolequal(min(p1),0, 1e-4);
mbg_asserttolequal(min(p2),0, 1e-4);