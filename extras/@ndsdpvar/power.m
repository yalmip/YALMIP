function X = power(X,d)
% UMINUS (overloaded)

s = X.dim;
X = sdpvar(X);
d = d(:);
X = power(X,d);
X = reshape(X,s);