function v = depends(X)

v = depends(X.cx);
for i = 1:length(X.P)
    v = [v depends(X.P{i})];
end
v = unique(v);
