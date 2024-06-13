function [V,D,permutation,failure] = dmpermblockeig(X,switchtosparse)
    
[permutation,aux1,aux2,blocks] = dmperm(X+speye(length(X)));
Xpermuted = X(permutation,permutation);

V = [];
D = [];
V = zeros(size(X,1),1);
top = 1;
left = 1;
anycholfail = 0;
failure = 0;

for i = 1:length(blocks)-1
    Xi = Xpermuted(blocks(i):blocks(i+1)-1,blocks(i):blocks(i+1)-1);
    [R,fail] = chol(Xi);
    anycholfail = anycholfail | fail;
    if fail && nnz(Xi)>0
        if length(Xi) >= switchtosparse           
            [vi,di,eigfail] = eigs(Xi,5,'SA');            
            if eigfail || isempty(di)
                res = 0;
                for j = 1:size(vi,2)
                    res(j) = norm(Xi*vi(:,j)-vi(:,j)*di(j,j));
                end
                % We only trust these
                notfailed = abs(res) <= 1e-12;
                vi = vi(:,notfailed);
                di = di(notfailed,notfailed);
                if length(vi) == 0
                    [vi,di,eigfail] = eigs(sparse(Xi),25,'SA');                    
                    if eigfail
                        res = 0;
                        for j = 1:size(vi,2)
                            res(j) = norm(Xi*vi(:,j)-vi(:,j)*di(j,j));
                        end
                        % We only trust these
                        notfailed = abs(res) <= 1e-12;
                        vi = vi(:,notfailed);
                        di = di(notfailed,notfailed);
                    end
                end
            end
        else
            [vi,di] = eig(full(Xi));
        end
        for j = 1:length(di)
            if di(j,j)<=0
                V(top:top+length(Xi)-1,left)=vi(:,j);
                left = left + 1;
                D = blkdiag(D,di(j,j));
            end
        end
    end
    top = top + length(Xi);
end


if (anycholfail && isempty(V)) || (anycholfail && all(diag(D)>0)) 
    % OK, we have a problem. The Cholesky factorization failed for some of
    % the matrices, but yet no eigenvalue decomposition revealed a negative
    % eigenvalue (due to convergence issues in the sparse eigs)
    failure = 1;
end