function Y = max(X)
% abs (overloaded)

narginchk(1,1);

d = size(X);
dnew = d;
dnew(1) = 1;
Y = reshape(max(reshape(X(:),d(1),[])),dnew);

