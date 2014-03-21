function test_cutsdp
c=[1 -2 3]';
x = binvar(3,1);
sdpvar t
A = [-1 2 0;-3 -4 1;0 0 -2];
P = sdpvar(3,3);
a = sdpvar(3,1);
F = [P >= 0, P*A'+A*P <= 0];
F = [F, a<=9*x];
F = [F, a>=-9*x];
F = [F, trace(P) == 1];
F = [F, [t (c+a)';c+a P] >=0];
ops = sdpsettings('solver','cutsdp','cutsdp.maxiter',10);
sol = solvesdp(F,t,ops)
% Just assert that it ran
mbg_asserttrue(sol.problem ~= -1);

binvar  x1 x2 x3 x4 x5 x6 x7 x8 x9 x10 x11 x12 x13 x14 x15
x = [x1 x2 x3 x4 x5 x6 x7 x8 x9 x10 x11 x12 x13 x14 x15]';
c = [-192 -96 -48 96 192 32 16 8 -16 -32 -80 -40 -20 40 80]';
A1 = -40*x1-20*x2-10*x3+20*x4+40*x5;
A2 = -16*x1-8*x2-4*x3+8*x4+16*x5+32*x11+16*x12+8*x13-16*x14-32*x15;
A3 = -23+32*x6+16*x7+8*x8-16*x9-32*x10+8*x11+4*x12+2*x13-4*x14-8*x15;
A = [A1 A2;A2 A3]
Cuts = diag(A) >= 0;
solvesdp(A>=0, c'*x,sdpsettings('solver','cutsdp'))
mbg_asserttrue(double(c'*x)==24);

solvesdp([A>=0, sum(x)==11], c'*x,sdpsettings('solver','cutsdp'))
mbg_asserttrue(double(c'*x)==24);

solvesdp([A>=0, A>= 0, sum(x)>=13], c'*x,sdpsettings('solver','cutsdp'))
mbg_asserttrue(double(c'*x)==68);


