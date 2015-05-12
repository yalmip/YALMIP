function issym=ishermitian(X)
%ISHERMITIAN Check if variable is Hermitian

n = X.dim(1);
m = X.dim(2);
issym = 0;
if (n==m)

    if isequal(X.conicinfo,[1 0])
        issym = 1;
        return
    elseif isequal(X.conicinfo,[-1 0])
        issym = 0;
    end
    
    if isa(X.basis,'lazybasis')
        issym = 1;
    else

        % What are the linear indicies to the transposed matrices
        %ind = reshape(reshape(1:n^2,n,n)',n^2,1);       
        [i,j,k] = find(X.basis);
        col = ceil(i/X.dim(1));
        row = i - (col-1)*X.dim(1);
        
        % Original linear indicies
        % row + (col-1)*X.dim(1)
        % Transposed linear indicies
        inew = col + (row-1)*X.dim(1);
        tranposedBasis = sparse(inew,j,k,size(X.basis,1),size(X.basis,2));
         
        if ~isreal(X.basis)
            residual = mid(X.basis-conj(tranposedBasis));
            %residual = mid(X.basis-conj(X.basis(ind,:)));
        else        
            residual = mid(X.basis-tranposedBasis);
            % Optimized for speed on huge sparse cases
            %[i,j,k] = find(X.basis);
            %basisTransposed = sparse(ind(i),j,k);
            %residual = mid(X.basis-basisTransposed);
            %%residual = mid(X.basis-X.basis(ind,:));
        end

        if nnz(residual)>0
            issym = norm(residual,1)<1e-10;
        else
            issym = 1;
        end
    end
end
