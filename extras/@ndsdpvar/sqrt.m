function X = abs(X)
% sqrt (overloaded)

X = reshape(sqrt(reshape(X,prod(X.dim),1)),X.dim);
