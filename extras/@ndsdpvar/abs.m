function X = abs(X)
%ABS (overloaded)

X = reshape(abs(reshape(X,prod(X.dim),1)),X.dim);
