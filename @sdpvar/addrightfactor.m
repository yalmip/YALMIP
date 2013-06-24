function X = addrightfactor(X,R)
if ~isempty(X.midfactors)
    if length(X.midfactors) >= 25
        % Just give up, this is not something which will be used in
        % STRUL...
        X = flush(X);
        return
    end
    for i = 1:length(X.midfactors)
        try
            X.rightfactors{i} = X.rightfactors{i}*R;
        catch
            X = flush(X);
            return
        end
    end
end
X = cleandoublefactors(X);
X = flushmidfactors(X);