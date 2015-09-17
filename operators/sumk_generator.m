function [Model,Properties] = sumk_generator(X,k,t)
[n,m] = size(X);
Z = sdpvar(n,m);
s = sdpvar(1,1);
if min(n,m)==1
    Model = (t-k*s-sum(Z) >= 0) + (Z >= 0) + (Z-X+s >= 0);
    Properties = struct('convexity','convex','monotonicity','increasing','definiteness','none','model','graph');
else
    Model = (t-k*s-trace(Z) >= 0) + (Z >= 0) + (Z-X+s*eye(n) >= 0);
    Properties = struct('convexity','convex','monotonicity','none','definiteness','none','model','graph');
end