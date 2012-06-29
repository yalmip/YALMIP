clear all;
close all;

N = 5;
n_lower = ones(N,1);
n_upper = 100*ones(N,1);
x_lower = ones(N,1);
x_upper = 4*ones(N,1);

h = [0 0 0 0 0 0]';
E = [1 0 0 0 0;
     0 1 1 0 0;
     0 0 1 0 0;
     0 1 1 1 0;
     0 0 1 1 0;
     0 0 0 1 0];

m = 2;
tl = 2;
Pl = 0.2;
NL = 50;
eta = 1.1;

DPx = [1000 1000 1000 1000 1000]';
d_lower = [1 1 1 1 1]';
d_upper = [1000 1000 1000 1000 1000]';
p0 = randn(5,1)+[0 0 0 0 0]';
d0 = 1+[0 0 0 0 0]';

Result = sdpvar(1,1);
Pt  = sdpvar(1,1);
n   = sdpvar(N,1);
x   = sdpvar(N,1);
D   = sdpvar(N,1);
P   = sdpvar(N,1);
fs  = sdpvar(1,1);

obj = (Result^m)*Pt;

F = set( Result^-1*(1+1*(h'*(E*n)))*fs^-1 <= 1 );
% for i = 1:N
%     F = F + set( n_lower(i) <= n(i) );
%     F = F + set( n(i) <= n_upper(i) );
%     F = F + set( x_lower(i) <= x(i) );
%     F = F + set( x(i) <= x_upper(i) );
%     F = F + set( n(i)^-1*D(i) + tl <= x(i)/fs );
%     F = F + set( (DPx(i)-p0(i)*d0(i))*D(i)^-1*P(i)^-1+p0(i)*P(i)^-1+d0(i)*D(i)^-1 <= 1 );
%     F = F + set( d_lower(i) <= D(i) );
%     F = F + set( D(i) <= d_upper(i) );
% end
    F = F + set( n_lower <= n );
    F = F + set( n <= n_upper );
    F = F + set( x_lower <= x );
    F = F + set( x <= x_upper );
    F = F + set( n.^-1.*D + tl <= x./fs );
    F = F + set( (DPx-p0.*d0).*(D.^-1).*(P.^-1)+p0.*(P.^-1)+d0.*(D.^-1) <= 1 );
    F = F + set( d_lower <= D);
    F = F + set( D <= d_upper);

% What I wanted to do in the previous version
F = F + set( Pt^-1*(x'*P)*fs + Pt^-1*((n.^eta)'*x)*NL*Pl*fs <= 1 );

% % Temporary solution
% temp = sdpvar(N,1);
% for i = 1:N
%     temp(i) = n(i)^eta;
% end    
% F = F + set( Pt^-1*(x'*P)*fs + Pt^-1*(temp'*x)*NL*Pl*fs <= 1 );

solvesdp(F,obj,sdpsettings('verbose',1))

disp('Result')
double(Result)
