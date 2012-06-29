function X = transposefactor(X)
if ~isempty(X.midfactors)
    for i = 1:length(X.midfactors)
        temp = X.leftfactors{i};
        X.leftfactors{i} = X.rightfactors{i}.';
        X.rightfactors{i} = temp.';
        X.midfactors{i} = X.midfactors{i}.';
    end
end
