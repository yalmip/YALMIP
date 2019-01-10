function X=permute(X,p)
%PERMUTE (overloaded)

if length(X.dim) < length(p)
    X.dim = [X.dim ones(1,length(p)-length(X.dim))];
end
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