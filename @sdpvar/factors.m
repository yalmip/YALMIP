function [L,M,R] = factors(X)

L = X.leftfactors;
M = X.midfactors;
R = X.rightfactors;
if isempty(M)
    return
end
% Prune. 
keep = 1:length(L);
for i = 1:length(L)
    if isnumeric(M{i})
        if nnz(M{i})==0
            keep(i) = 0;
        end
    end
end
keep = find(keep);
L = {L{keep}};
M = {M{keep}};
R = {R{keep}};
