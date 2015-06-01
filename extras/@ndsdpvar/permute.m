function X=permute(X,p)
%PERMUTE (overloaded)

i = 1:prod(X.dim);
i = reshape(i,X.dim);
i = permute(i,p);
X.basis = X.basis(i,:);
X.dim = X.dim(p);
while length(size(X)) > 2 && X.dim(end) == 1;
    X.dim = X.dim(1:end-1);
    if length(X.dim)==2
        X = reshape(sdpvar(X),X.dim);
    end
end