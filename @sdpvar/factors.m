function [L,M,R] = factors(X)

% Author Johan Löfberg
% $Id: factors.m,v 1.2 2009-10-14 09:10:35 joloef Exp $

L = X.leftfactors;
M = X.midfactors;
R = X.rightfactors;
if isempty(M)
    return
end
% Prune. 
keep = 1:length(L);
for i = 1:length(L)
    if isa(M{i},'double')
        if nnz(M{i})==0
            keep(i) = 0;
        end
    end
end
keep = find(keep);
L = {L{keep}};
M = {M{keep}};
R = {R{keep}};
