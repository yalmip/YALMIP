function X = sum(varargin)
% SUM (overloaded)

Y = varargin{1};
X = Y;
X.basis = [];
if nargin == 2 && isequal(varargin{2},length(X.dim))
    % smash slices
    X.basis = kron(ones(1,X.dim(end)),speye(prod(X.dim(1:end-1))))*Y.basis;
    X.dim = X.dim(1:end-1);
else
    for i = 1:size(Y.basis,2)
        base = reshape(full(Y.basis(:,i)),X.dim);
        base = sum(base,varargin{2:end});
        X.basis = [X.basis sparse(base(:))];
    end
    X.dim = size(base);
end
X.conicinfo = [0 0];
X = clean(X);