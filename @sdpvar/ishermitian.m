function issym=ishermitian(X)
%ISHERMITIAN Check if variable is Hermitian

% Author Johan Löfberg
% $Id: ishermitian.m,v 1.8 2006-12-13 13:13:32 joloef Exp $

n = X.dim(1);
m = X.dim(2);
issym = 0;
if (n==m)

    if isequal(X.conicinfo,[1 0])
        issym = 1;
        return
    end
    
    if isa(X.basis,'lazybasis')
        issym = 1;
    else

        % What are the linear indicies to the transposed matrices
        ind = reshape(reshape(1:n^2,n,n)',n^2,1);

        if ~isreal(X.basis)
            residual = mid(X.basis-conj(X.basis(ind,:)));
        else
            residual = mid(X.basis-X.basis(ind,:));
        end

        if nnz(residual)>0
            issym = norm(residual,1)<1e-10;
        else
            issym = 1;
        end
    end
end
