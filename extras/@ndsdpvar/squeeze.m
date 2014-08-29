function X = squeeze(X)
% squeeze (overloaded)

dummy = reshape(ones(prod(X.dim),1),X.dim);
dummy = squeeze(dummy);
X.dim = size(dummy);
X = clean(X);
