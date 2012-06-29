function X=transpose(X)
%TRANSPOSE (overloaded)

% Author Johan Löfberg 
% $Id: transpose.m,v 1.1 2006-08-10 18:00:23 joloef Exp $

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

