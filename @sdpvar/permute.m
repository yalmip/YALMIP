function X=permute(X,p)
%PERMUTE (overloaded)

% Author Johan Löfberg 
% $Id: permute.m,v 1.1 2007-02-14 16:12:20 joloef Exp $   

i = 1:prod(X.dim);
i = reshape(i,X.dim);
i = permute(i,p);
X.basis = X.basis(i,:);
X.dim = X.dim(p);