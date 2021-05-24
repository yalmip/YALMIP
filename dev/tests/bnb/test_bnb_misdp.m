function tests = test_bnb_misdp
tests = functiontests(localfunctions);

function test_misdp(testCase)
% Test two regression bugs
% 1. quadratic costs in binary SDP
% 2. OR on SDP constraints

X = sdpvar(3,3);
obj = trace((X-2*eye(3))*(X-2*eye(3))');
sol = optimize((X<=3*eye(3)) + ((X>=eye(3)) | (X<=-eye(3))) + (-50 <= X(:) <= 50),obj,sdpsettings('verbose',0,'debug',1,'sedumi.free',0));

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
