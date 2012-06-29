x1 = sdpvar(1,1);
x2 = sdpvar(1,1);
x3 = sdpvar(1,1);
x4 = sdpvar(1,1);
x5 = sdpvar(1,1);
x6 = sdpvar(1,1);
x7 = sdpvar(1,1);
x8 = sdpvar(1,1);
x9 = sdpvar(1,1);
x10 = sdpvar(1,1);
x11 = sdpvar(1,1);
x12 = sdpvar(1,1);
x13 = sdpvar(1,1);
x14 = sdpvar(1,1);
x15 = sdpvar(1,1);

x = [x1;x2;x3;x4;x5;x6;x7;x8;x9;x10;x11;x12;x13;x14;x15];

c = [-192 -96 -48 96 192 32 16 8 -16 -32 -80 -40 -20 40 80]';
A1 = -40*x1-20*x2-10*x3+20*x4+40*x5;
A2 = -16*x1-8*x2-4*x3+8*x4+16*x5+32*x11+16*x12+8*x13-16*x14-32*x15;
A3 = -23+32*x6+16*x7+8*x8-16*x9-32*x10+8*x11+4*x12+2*x13-4*x14-8*x15;

A = [A1 A2;A2 A3];
F = lmi(A>0) + lmi(x>0)+lmi(x<1)

solvesdp(F+binary(x),c'*x,sdpsettings('solver','cp'));


F_lin = set(diag(A))+lmi(x>0)+lmi(x<1);


solvesdp(F_lin,c'*x,sdpsettings('solver','glpk'))

for i = 1:5
    [v,d] = eig(double(A));
    F_lin = F_lin + set(v(:,1)'*A*v(:,1) > 0);
    solvesdp(F_lin+integer(x),c'*x,sdpsettings('solver','glpk'));
    double(c'*x)
end



solvesdp(F,[],c'*x,sdpsettings('verbose',2,'bnb.maxiter',25))


double(c'*x)
double(x)
