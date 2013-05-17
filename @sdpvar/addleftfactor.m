function X = addleftfactor(X,L)
if ~isempty(X.midfactors)
    for i = 1:length(X.midfactors)
        try
            X.leftfactors{i} = L*X.leftfactors{i};
        catch
            X = flush(X);
            return
        end
    end
end
X = cleandoublefactors(X);
X = flushmidfactors(X);
