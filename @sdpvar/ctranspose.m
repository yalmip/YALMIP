function X=ctranspose(X)
%CTRANSPOSE (overloaded)

% Author Johan Löfberg
% $Id: ctranspose.m,v 1.7 2009-08-26 01:46:41 joloef Exp $

if isa(X,'blkvar')
    X = sdpvar(X);
end

n = X.dim(1);
m = X.dim(2);
ind = reshape(reshape(1:n*m,n,m)',n*m,1);
if isreal(X.basis)
    X.basis = X.basis(ind,:);
else
    X.basis = conj(X.basis(ind,:));
end
X.dim(1) = m;
X.dim(2) = n;
% Reset info about conic terms
X.conicinfo = [0 0];
X = ctransposefactor(X);