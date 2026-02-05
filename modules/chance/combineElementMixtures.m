function distribution = combineElementMixtures(distribution,id)

if isempty(distribution.mixture)
    return
end

% Draw a mixture for every unique component (the problem is multivariate
% mixtures)
[unique_distributionIDs,index_to,index_from] = unique(id,'stable') ;

Combs = allSequences(length(unique_distributionIDs),size(distribution.mixture,2));


for k = 2:length(distribution.parameters)
    newParams = {};  
    newParam = distribution.parameters{k}{1};
    for i = 1:length(Combs)
          
    end
end

return


function C = allSequences(m,n)
% Returns an (m^n) x n matrix whose rows are all length-n sequences from 1:m.
vectors = cell(1,n);
[vectors{:}] = ndgrid(1:m);
C = cell2mat( cellfun(@(V) V(:), vectors, 'UniformOutput', false) );

