function pass = akerblad

G=tf(1,[1 1 1]);
Gd=c2d(G,0.1);
[A,B]=ssdata(Gd);  
[n,m]=size(B);
N=100;

Ci=[1 0;0 1; 0 0];
Di=[0; 0; 1];
p=3;
q=4;
Ei=[1 3;-1 -3;0 0; 0 0];
Fi=[0;0;1; -1];
bi=[10;10;1;1];
P=[1 0;0 1;];
W=[1 0;0 1;];
E=[];F=[];C=[];D=[];
ci=[1;1;zeros(p,1)];
x0=[2;2];
b=[];
c=[];
for k=1:N,
  E(:,:,k)=Ei;
  F(:,:,k)=Fi;
  C(:,:,k)=Ci;
  D(:,:,k)=Di;
  b=[b;bi];
  c=[c;ci];
end
c=[c;1;1;zeros(n,1);1;zeros(n,1)];
d=sparse([],[],[],(N+1)*n,1);
d(1:n)=x0;
H=sparse([],[],[]);
H(1:n,1+1:n+1)=eye(n);
M=sparse([],[],[],N*q,N*(n+m+1)+1+n);
L=sparse([],[],[]);

for k=1:N,
  M(k*q+1-q:k*q,k*(n+m+1)-n-m:k*(n+m+1))=[zeros(q,1) E(:,:,k) F(:,:,k)];
  L(k*(p+2)+1-p-2:k*(p+2),k*(n+m+1)-n-m:k*(n+m+1))=[-1 zeros(1,n+m);1 zeros(1,n+m);zeros(p,1) -2*C(:,:,k) -2*D(:,:,k)];
  H(k*n+1:k*n+n,k*(n+m+1)+1-n-m:k*(1+n+m)+n+1)=[-A -B zeros(n,1) ...
		    eye(n)];
end
L(k*(p+2)+1:(k)*(p+2)+2+n+1+n,k*(n+m+1)+1:k*(n+m+1)+n+1)=[-1 ...
		    zeros(1,n);1 zeros(1,n);zeros(n,1) -2*P;zeros(1,n+ ...
						  1);zeros(n,1) -W];
f=sparse([],[],[],N*(n+m+1)+n+1,1);
f(1:n+m+1:N*(n+m+1)+n+1)=1;
[yy,ss,rr,sol]=my_socp2(f,d,c,b,L,M,H,N);

norm(yy)
norm(ss)









function [y,s,r,sol]=my_socp2(f,b,c,d,A,B,H,N)
%[t,x,s,r,gamma,lambda]=my_socp(f,c,d,A,B)
%Solves the primal and dual LP+SOCP on the form
%min f'y                |  max  -d'gamma-c'lambda-b'rho
%st  s=c-Ay             |  st   -A'lambda-B'gamma-H'rho=f
%    r=d-By             |       0<gamma, ||lambda_1||<lambda_0 
%    Hy=b               |  
%    0<r, ||s_1||<s_0   |
[n,m]=size(A);
[p,m]=size(B);
[q,m]=size(H);

tic
yalmip('clear')

y=sdpvar(m,1);
s = c-A*y;
r = d-B*y;
F=lmi('0<r');
F=addlmi(F,'H*y=b');
for i=1:N
F=addlmi(F,'||s(i*(2+3)-3:i*5)||<s(i*5-4)');
end
F=addlmi(F,'||s((i+1)*(2+3)-3:(i+1)*5-1)||<s((i+1)*5-4)');
F=addlmi(F,'||s(n-1:n)||<s(n-2)');
toc
sol=solvesdp(F,[],f'*y,sdpsettings('sedumi.beta',0.8))

y=double(y); s=double(s); r=double(r);



