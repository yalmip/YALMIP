function [v,Model] = eigv(X)
% v(i) epigraph of i-th largest eigenvalue of X
[n,m] = size(X);
v = sdpvar(n,1);
Model = [trace(X)==sum(v), v(1:end-1) >= v(2:end)];
for k = 1:n-1             
    Z = sdpvar(n);
    s = sdpvar(1);
    t = sum(v(1:k));
    Model = [Model, (t-k*s-trace(Z) >= 0) + (Z >= 0) + (Z-X+s*eye(n) >= 0)];    
end