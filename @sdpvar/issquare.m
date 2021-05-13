function issquare=issquare(X)
%ISSQUARE Check if variable is square

n = X.dim(1);
m = X.dim(2);
issquare = (n==m);
