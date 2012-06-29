function X=transpose(X)
%TRANSPOSE (overloaded)

% Author Johan Löfberg 
% $Id: transpose.m,v 1.7 2006-07-26 20:17:58 joloef Exp $

if isa(X,'blkvar')
    X = sdpvar(X);
end

n = X.dim(1);
m = X.dim(2);
ind = reshape(reshape(1:n*m,n,m)',n*m,1);
X.basis = X.basis(ind,:);
X.dim(1) = m;
X.dim(2) = n;
% Reset info about conic terms
X.conicinfo = [0 0];
X = transposefactor(X);

