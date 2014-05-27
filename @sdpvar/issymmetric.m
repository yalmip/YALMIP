function issym=issymmetric(X)
%ISSYMMETRIC Check if variable is symmetric

n = X.dim(1);
m = X.dim(2);
if (n==m)
    % What are the linar indicies to the transposed matrices
    if isa(X.basis,'lazybasis')
        issym = 1;
    else
        ind = reshape(reshape(1:n^2,n,n)',n^2,1);
        % We use mid in case there is interval data
        issym = norm(mid(X.basis-X.basis(ind,:)),1)<1e-10;   
    end
else
    issym = 0;
end

