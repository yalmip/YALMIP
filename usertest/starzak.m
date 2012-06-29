function fail = starzak

N = 15;  % Number of balls
n = 6;   % ...in R^n

% Generate som balls
randn('seed',123456);
for i = 1:N
  P = randn(n);P = P*P'+eye(n)*0.1;
  Pi{i}=P;
end

% Define LMIs from S-procedure
tau = sdpvar(N,1);
P = sdpvar(n,n);
G = lmi('P>0');
F = lmi(P > 1e3*eye(n)*sqrt(eps)); % Numericl reasons
for i = 1:N
  F = F+lmi(tau(i)*Pi{i} > P);
end
F = F+lmi(tau<1);

% MAXDET
sol = solvesdp(F,-logdet(P),sdpsettings('verbose',0));
fail=getfail(sol.problem,-log(det(double(P))),10.8514,checkset(F));