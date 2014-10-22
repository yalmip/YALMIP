function [V,D,permutation] = blockeig(X,switchtosparse)
    
[permutation,aux1,aux2,blocks] = dmperm(X);
Xpermuted = X(permutation,permutation);

V = [];
D = [];
V = zeros(size(X,1),1);
top = 1;
left = 1;

for i = 1:length(blocks)-1
    Xi = Xpermuted(blocks(i):blocks(i+1)-1,blocks(i):blocks(i+1)-1);
    [R,fail] = chol(Xi);
    if fail
        if length(Xi) >= switchtosparse                                            
            [vi,di] = eigs(sparse(Xi),5,'SA');
        else
            [vi,di] = eig(Xi);
        end
        for j = 1:length(di)
            if di(j,j)<0
                V(top:top+length(Xi)-1,left)=vi(:,j);
                left = left + 1;
                D = blkdiag(D,di(1,1));
            end
        end
    end
    top = top + length(Xi);
end
