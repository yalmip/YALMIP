function X = flushmidfactors(X)
for i = 1:length(X.midfactors)
    if isa(X.midfactors{i},'sdpvar')
        X.midfactors{i} = flush(X.midfactors{i});
    end
end