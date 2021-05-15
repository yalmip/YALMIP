function tests = test_bnb
tests = functiontests(localfunctions);

function test_milp(testCase)
rand('seed',1234);

a = [1 2 3 4 5 6]';
t = (0:0.02:2*pi)';
x = [sin(t) sin(2*t) sin(3*t) sin(4*t) sin(5*t) sin(6*t)];
y = x*a+(-4+8*rand(length(x),1));

a_hat = intvar(6,1);

residuals = y-x*a_hat;
bound = sdpvar(length(residuals),1);
F = (-bound <= residuals <= bound);
ops = sdpsettings('solver','bnb','verbose',0);

% Test QP
obj = sum(bound);
sol = optimize(F,obj,ops);
testCase.assertTrue(sol.problem == 0);
testCase.assertTrue(abs(value(obj) - 6.168422746718130e+002) <= 1e-5);


function test_minlp(testCase)
randn('seed',12345);
rand('seed',12345);

x = sdpvar(5,1);
A = randn(15,5);
b = rand(15,1)*10;

obj = sum(x) + sum((x-3).^4);
constraints = (A*x <= b) + (integer(x));
sol = optimize(constraints,obj,sdpsettings('solver','bnb','bnb.solver','fmincon','verbose',0));

testCase.assertTrue(sol.problem == 0);


function test_miqp(testCase)
rand('seed',1234);

a = [1 2 3 4 5 6]';
t = (0:0.02:2*pi)';
x = [sin(t) sin(2*t) sin(3*t) sin(4*t) sin(5*t) sin(6*t)];
y = x*a+(-4+8*rand(length(x),1));

a_hat = intvar(6,1);

residuals = y-x*a_hat;
bound = sdpvar(length(residuals),1);
F = (-bound <= residuals <= bound);
ops = sdpsettings('solver','bnb','verbose',0);

% Test QP
obj = residuals'*residuals;
sol = optimize((residuals <= 50),obj,ops);
testCase.assertTrue(sol.problem == 0);
testCase.assertTrue(abs(value(obj) - 1.605058709613011e+003) <= 1e-5);


function test_misdp(testCase)
% Test two regression bugs
% 1. quadratic costs in binary SDP
% 2. OR on SDP constraints

X = sdpvar(3,3);
obj = trace((X-2*eye(3))*(X-2*eye(3))');
sol = optimize((X<=3*eye(3)) + ((X>=eye(3)) | (X<=-eye(3))) + (-50 <= X(:) <= 50),obj,sdpsettings('verbose',0));

testCase.assertTrue(sol.problem==0);
testCase.assertTrue(abs(value(obj)) <= 1e-5);


Loc_n=4;
Sen_n=3;
Nvar=Loc_n*Sen_n;  %Variable number
%Constants
C=[1 1;1 0;0 1;0 1];
s1=5.0693;
s4=3.8259;
Cov_M=[1.501 0.523 0.978 0.978; 3.002 1.046 1.956 1.956; 4.503 1.569 2.934 2.934]';
Z=[2500 1500 800]; %$,cost
%variable claim
X=binvar(Nvar,1);    % [q11 q12 q13 q21 q22 q23 q31 q32 q33 q41 q42 q43]
%Constraints (LMIs)
for i=1:Loc_n
    k=3*i-2;
    Xigma_inv(i,i)=1./Cov_M(i,:).^2*X(k:k+2);
end
P1(2:3,2:3)=C'*Xigma_inv*C;
P1(1,1)=s1;
P1(1,2:3)=C(1,:);
P1(2:3,1)=C(1,:)';
P2=P1;
P2(1,1)=s4;
P2(1,2:3)=C(4,:);
P2(2:3,1)=C(4,:)';
P3(1,1)=1-sum(X(1:3));
P3(2,2)=1-sum(X(4:6));
P3(3,3)=1-sum(X(7:9));
P3(4,4)=1-sum(X(10:12));
F=(P1>=0)+(P2>=0)+(diag(P3)>=0);
cc=zeros(1,Nvar);
cc(1:Sen_n)=Z;
cc(Sen_n+1:2*Sen_n)=Z;
cc(2*Sen_n+1:3*Sen_n)=Z;
cc(3*Sen_n+1:end)=Z;
obj=cc*X;
sol = optimize(F,obj,sdpsettings('solver','bnb','verbose',0));


function test_micones(testCase)
randn('seed',123456);
n = 5;
%  65.63177408018771   6.95052852960134 -58.66230788287043 -70.86445052912004  58.37564024957875
%  17.58996579270265  65.53001366712473  18.43985990767680 -58.43579533652219 -72.32340882851518
% -46.45728349558073  15.51300078651673  68.32913865172714  18.82218749426777 -57.90310836270277
% -79.19194296790369 -53.92234161527504  21.23787385808123  65.53955212353939  17.48548000466040
%  51.25486176425121 -73.14184210768217 -50.35188380380193  19.83022901147247  67.22698658022694
P = toeplitz(randn(n,1)*100)+randn(n,n)*5;
Z = intvar(n,n,'toeplitz');
t = sdpvar(n,n,'full');
e = P(:)-Z(:);
ops = sdpsettings('solver','bnb','verbose',0);

F = (-t(:) <= P(:)-Z(:) <= t(:));
obj = sum(sum(t));
sol = optimize(F,obj,ops);
testCase.assertTrue(sol.problem == 0);
testCase.assertTrue(abs(value(obj) - 66.18236738983525) <=  1e-4);

F = ([]);
obj = norm(e,1);
sol = optimize(F,obj,ops);
testCase.assertTrue(sol.problem == 0);
testCase.assertTrue(abs(value(obj) - 66.18236738983525) <= 1e-4);

obj = e'*e;
F = ([]);
sol = optimize(F,obj,ops);
testCase.assertTrue(sol.problem == 0);
testCase.assertTrue(abs(value(obj) - 3.352603490492911e+002) <= 1e-4);

t = sdpvar(1,1);
obj = t;
F = (cone(e,t));
sol = optimize(F,obj,ops);
testCase.assertTrue(sol.problem == 0);
testCase.assertTrue(abs(value(obj) - 18.31011603130778) <= 1e-4);

t = sdpvar(1,1);
obj = norm(e);
F = ([]);
sol = optimize(F,obj,ops);
testCase.assertTrue(sol.problem == 0);
testCase.assertTrue(abs(value(obj) - 18.31011603130778) <= 1e-4);

obj = t;
F = ([t e';e eye(length(e))]>=0);
sol = optimize(F,obj,ops);
testCase.assertTrue(sol.problem == 0);
testCase.assertTrue(abs(value(obj) - 3.352603420494530e+002) <= 2e-4);



function test_migp(testCase)
x = sdpvar(7,1);

% Data
a     = ones(7,1);
alpha = ones(7,1);
beta  = ones(7,1);
gamma = ones(7,1);
f = [1 0.8 1 0.7 0.7 0.5 0.5]';
e = [1 2 1 1.5 1.5 1 2]';
Cout6 = 10;
Cout7 = 10;

% Model
C = alpha+beta.*x;
A = sum(a.*x);
P = sum(f.*e.*x);
R = gamma./x;

D1 = R(1)*(C(4));
D2 = R(2)*(C(4)+C(5));
D3 = R(3)*(C(5)+C(7));
D4 = R(4)*(C(6)+C(7));
D5 = R(5)*(C(7));
D6 = R(6)*Cout6;
D7 = R(7)*Cout7;

% Constraints
F = (x >= 1) + (P <= 20) + (A <= 100);

% Objective
D = max([(D1+D4+D6),(D1+D4+D7),(D2+D4+D6),(D2+D4+D7),(D2+D5+D7),(D3+D5+D6),(D3+D7)]);

% Solve integer problem
ops = sdpsettings('solver','bnb','verbose',0);
sol = optimize(F+(integer(x)),D,ops);

testCase.assertTrue(sol.problem == 0);
testCase.assertTrue(all(abs(value(x) - [ 2     3     3     3     2     3     3]') <= 1e-3));
testCase.assertTrue(abs(value(D)-(8+1/3)) <= 1e-3);


function test_diw(testCase)
[F,h] = loadsdpafile('diw_15.dat-s');
sol = optimize(F,h,sdpsettings('solver','bnb'));
testCase.assertTrue(sol.problem == 0);
testCase.assertTrue(abs(value(h) - -95) <= 1e-2);
function test_colon(testCase)
[F,h] = loadsdpafile('coloncancer_1_100_5.dat-s')
sol = optimize(F,h,sdpsettings('solver','bnb'));
testCase.assertTrue(sol.problem == 0);
testCase.assertTrue(abs(value(h) - 127.47) <= 1e-2);
function test_clique(testCase)
[F,h] = loadsdpafile('clique_20_k3_6_7.dat-s')
sol = optimize(F,h,sdpsettings('solver','bnb'));
testCase.assertTrue(sol.problem == 0);
testCase.assertTrue(abs(value(h) - 147) <= 1e-2);
function test_bar(testCase)
[F,h] = loadsdpafile('2x3_3bars.dat-s')
sol = optimize(F,h,sdpsettings('solver','bnb'));
testCase.assertTrue(sol.problem == 0);
testCase.assertTrue(abs(value(h) - 2.118) <= 1e-2);