function Q=getbasematrix(X,ind)
%GETBASEMATRIX Internal function to extract basematrix for variable IND

if ind(1) == 0
    base = X.basis(:,1);
    Q = reshape(base,X.dim(1),X.dim(2));
else
    here = find(X.lmi_variables==ind(1));
    if isempty(here)
        Q = sparse(X.dim(1),X.dim(2));
    else
        base = X.basis(:,here+1);
        Q = reshape(base,X.dim(1),X.dim(2));
    end
end
if length(ind)>1
    Q = [Q getbasematrix(X,ind(2:end))];
end


