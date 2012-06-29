function [H,S]=create_CHS(A,B,C,N,M);
% CREATE_CHS Internal function to define matrices for MPC problem

% Author Johan Löfberg
% $Id: create_CHS.m,v 1.2 2004-07-02 08:17:29 johanl Exp $

% How many states
n=length(A);p=size(C,1);
[dummy,m]=size(B);

H=zeros(N*p,length(A));
S=zeros(N*p,N*m);

Acum=A;
for j=1:N,
  H(1+p*(j-1):p*j,:)=C*Acum;
  Acum=Acum*A;
end;

for j=1:N,
  Acum=eye(n);
  for k=j:-1:1,
    S(1+p*(j-1):p*j,1+m*(k-1):k*m)=C*Acum*B;
    Acum=Acum*A;
  end;
end;
