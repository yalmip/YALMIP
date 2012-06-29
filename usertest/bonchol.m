clear all

yalmip('clear')
alfa=1;
n = 6;

load Q.mat
w  = -0.5;
Qg = frsp(Qde,10^w);
M0 = Qg(1:length(Qg)-1,1:length(Qg)-1);

d1 = sdpvar(1,1); d2 = sdpvar(1,1); d3 = sdpvar(1,1); d4 = sdpvar(1,1); d5 = sdpvar(1,1);
d = diag([d1 d2 d3 d4 d5 d5]); 

h1 = sdpvar(1,1); h2 = sdpvar(1,1); h3 = sdpvar(1,1); h4 = sdpvar(1,1); 
h = diag([h1 h2 h3 h4 0 0]); 

o = zeros(2*n,2*n);
I = eye(2*n);

M1 = real(M0); M2 = imag(M0);
M  = [M1 M2; -M2 M1];
jM = [-M2 M1;-M1 -M2];

% Define D & increase size to account for complex numbers
D1 = d; D2 = zeros(n,n);
D = [D1 D2; -D2 D1];

% Define H such that G=D'H and increase size to account for complex numbers
H1 = h; H2 = zeros(n,n); 
H = [H1 H2; -H2 H1];
jH = [-H2 H1;-H1 -H2];

W = [ I         D*M-jH     o        o;
     M'*D'+jH'    o      alfa*D'    H';
      o         alfa*D    -I        o;
      o           H        o       -I];
    
Y = sdpvar(8*n,8*n,'symmetric');
Z = sdpvar(8*n,8*n,'symmetric');

L = lmi(d-1e-3);
L = addlmi(L,[Y W;W' Z]);

objective = trace(Y)+trace(Z);

options = sdpsettings('Solver','sedumi','sedumi.maxiter',2);
options.sedumi.alg   = 2;
options.sedumi.theta = 0.5;
options.sedumi.beta  = 0.5;
options.sedumi.eps   = 1e-3;

solution = solvesdp(L,[],objective,options);

W = double(W);
D = double(D);
D = D(1:n,1:n);
H = double(H);
H = H(1:n,1:n);
M = M1+j*M2;
G = D'*H;
    
MM = M'*D*M+j*(G*M-M'*G)-alfa*D;

lambda = eig(MM);