function X = addleftfactor(X,L)
if ~isempty(X.midfactors)
    if length(X.midfactors) >= 25
        % Just give up, this is not something which will be used in
        % STRUL...
        X = flush(X);
        return
    end
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
