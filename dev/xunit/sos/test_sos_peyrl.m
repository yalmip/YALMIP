function tests = test_sos_peyrl
tests = functiontests(localfunctions);

function test1(dummy)

sdpvar x1 x2 x3 x4;
x = [x1;x2;x3;x4];
n = length(x);

% System matrices:
Ao = [0  1  0 0;
     0  0  1 0;
     0  0  0 1;
     -1 0 -2 0];

bo = [0 0 0 1]';

fo = [-1 -1 -2 -1];

% Transformation-Matrix:
Pinv = [1  -1 -1  0;
        -1  1 -1  0;
        0  -1  1  0;
        0   2  0 -1];

A = Ao;%inv(Pinv)*Ao*Pinv;
b = bo;%inv(Pinv)*bo;

f = fo;%fo*Pinv;

% System equations:
f1 = A*x+b*f*x;
f2 = A*x+b;
f3 = A*x-b;

% Cell boundaries:
g11 = f*x+1;
g12 = -f*x+1;
g21 = f*x-1;
g31 = -f*x-1;
h120 = f*x-1;
h130 = f*x+1;


% The Lyapunov functions Vi(x):
disp('Constructing Lyapunov functions')

N = 4;

% The Lyapunov functions Vi(x):
lm = monolist(x,N,2);
p1 = sdpvar(size(lm,1),1); define(p1);
V1 = p1'*lm;

lm = monolist(x,N);
p2 = sdpvar(size(lm,1),1); define(p2);
V2 = p2'*lm;
p3 = sdpvar(size(lm,1),1); define(p3);
V3 = p3'*lm;


sdpvar eps1 eps2 eps3 eps4 eps5 eps6
pd1 = eps1*sum(x.^N);
pd2 = eps2*sum(x.^N);
pd3 = eps3*sum(x.^N);
pd4 = eps4*sum(x.^N);
pd5 = eps5*sum(x.^N);
pd6 = eps6*sum(x.^N);

% Create the aijk, bijk SOS-polys and cij:
lm = monolist(x,N/2-1);
A11 = sdpvar(size(lm,1)); (A11); a11 = lm'*A11*lm;
A12 = sdpvar(size(lm,1)); (A12); a12 = lm'*A12*lm;
A21 = sdpvar(size(lm,1)); (A21); a21 = lm'*A21*lm;
A31 = sdpvar(size(lm,1)); (A31); a31 = lm'*A31*lm;

B11 = sdpvar(size(lm,1)); (B11); b11 = lm'*B11*lm;
B12 = sdpvar(size(lm,1)); (B12); b12 = lm'*B12*lm;
B21 = sdpvar(size(lm,1)); (B21); b21 = lm'*B21*lm;
B31 = sdpvar(size(lm,1)); (B31); b31 = lm'*B31*lm;

lm = monolist(x,N-1);
C12 = sdpvar(size(lm,1)); (C12); c12 = lm'*C12*lm;
C13 = sdpvar(size(lm,1)); (C13); c13 = lm'*C13*lm;

% Constraints:
F = (sos(V1 - a11*g11 - a12*g12));
F = F + (sos(V2 - a21*g21 ));
F = F + (sos(V3 - a31*g31 ));

F = F + (sos(-jacobian(V1,x)*f1 - b11*g11 - b12*g12 ));
F = F + (sos(-jacobian(V2,x)*f2 - b21*g21) );
F = F + (sos(-jacobian(V3,x)*f3 - b31*g31) );

F = F + (coefficients(V1+c12*h120-V2,x)==0);
F = F + (coefficients(V1+c13*h130-V3,x)==0);


F = F + (eps1>=0.1) + (eps2>=0.1) + (eps3>=0.1);
F = F + (eps4>=0.1) + (eps5>=0.1) + (eps6>=0.1);

F = F + (A11>=0) + (A12>=0) + (A21>=0) + (A31>=0) + (B11>=0) + ...
    (B12>=0) + (B21>=0) + (B31>=0); 

% Call solver:
parametric=recover(setdiff(depends(F),depends(x)));
[sol,v,Q,residual,model]=solvesos(F,[],sdpsettings('sos.post',1,'sedumi.free',0),parametric);

assert(sol.problem==0 | sol.problem == 4);
assert(norm(residual) < 1e-5);
