function X = subsref(Y,refs)
% SUBSREF (overloaded)

% Author Johan Löfberg
% $Id: subsref.m,v 1.7 2006-07-13 20:37:30 joloef Exp $

X = Y;
base = reshape(1:size(Y.basis,1),X.dim);
base = subsref(base,refs);
X.basis = X.basis(base(:),:);
X.dim = size(base);
X = clean(X);