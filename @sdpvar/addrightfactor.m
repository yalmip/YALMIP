function X = addrightfactor(X,R)
if ~isempty(X.midfactors)
    for i = 1:length(X.midfactors)
        try
            X.rightfactors{i} = X.rightfactors{i}*R;
        catch
            X = flush(X);
        end
    end
end
X = cleandoublefactors(X);
X = flushmidfactors(X);