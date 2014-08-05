function test_dualize_sdp_11

yalmip('clear')
N = 2;
X = sdpvar(N,N,'hermitian','complex');
x = sdpvar(3,1);
t = sdpvar;
obj = t;
F = [X>=0,x>=1, X(1)+x(3) == 7];
F = [F,cone([1;t])];
opts = sdpsettings('dualize',0);
sol1 = solvesdp(F,-obj,opts);
o1 = double(obj);
opts = sdpsettings('dualize',1);
sol2 = solvesdp(F,-obj,opts);
o2 = double(obj);
mbg_asserttolequal(sol1.problem,0);
mbg_asserttolequal(sol2.problem,0);
mbg_asserttolequal(o1,o2,1e-4);
