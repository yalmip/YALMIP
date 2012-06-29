function X = abs(X)
% sqrt (overloaded)

% Author Johan Löfberg
X = reshape(sqrt(reshape(X,prod(X.dim),1)),X.dim);
