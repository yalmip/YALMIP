function Q=getbasematrix(X,ind)
%GETBASEMATRIX Internal function to extract basematrix for variable IND

if ind==0
    base = X.basis(:,1);
    Q = reshape(base,X.dim);
    return;
end

here = find(X.lmi_variables==ind);
if isempty(here)
    error
else
    base = X.basis(:,here+1);
    Q = reshape(base,X.dim);
end


