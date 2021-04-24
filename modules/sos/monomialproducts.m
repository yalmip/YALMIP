function [newton_m2,N_unique,newton_m2_unique] = monomialproducts(N,n);
%MONOMIALPRODUCTS  Internal function used for monomial reduction

% Exponents in squared monomials

N_unique = [];
for i = 1:size(N,1)
    newton_m2{i} = [];
    n = size(N{i},1);
    for j = 1:n
        newton_m2{i} = [newton_m2{i};[(1:n)' repmat(j,n,1) N{i}(1:n,:)+repmat(N{i}(j,:),n,1)]];
    end
    % Whoops, double copies of diagonal (we want double copies of non-diagonals though)
    if isempty(newton_m2{i})
        newton_m2_unique{i} = [];
    else
        [dummy,j,dummy2] = uniquesafe(newton_m2{i}(:,1:2),'rows');
        newton_m2{i} = newton_m2{i}(j,:);
        % Extract unique monomial products
        [dummy,j,dummy2] = uniquesafe(newton_m2{i}(:,3:end),'rows');
        newton_m2_unique{i} = newton_m2{i}(j,:);
    end
    N_unique = [N_unique;newton_m2_unique{i}];   
end
if ~isempty(N_unique)
    [dummy,j,dummy2] = uniquesafe(N_unique(:,3:end),'rows');
    N_unique = N_unique(j,:);
end
