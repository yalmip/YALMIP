function X=permute(X,p)
%PERMUTE (overloaded)

i = 1:prod(X.dim);
i = reshape(i,X.dim);
i = permute(i,p);
X.basis = X.basis(i,:);
X.dim = X.dim(p);