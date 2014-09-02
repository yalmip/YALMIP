function X = subsref(Y,refs)
% SUBSREF (overloaded)

X = Y;
base = reshape(1:size(Y.basis,1),X.dim);
base = subsref(base,refs);
X.basis = X.basis(base(:),:);
X.dim = size(base);
X = clean(X);