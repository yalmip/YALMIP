function X = abs(X)
% abs (overloaded)

% Author Johan Löfberg
X = reshape(abs(reshape(X,prod(X.dim),1)),X.dim);
