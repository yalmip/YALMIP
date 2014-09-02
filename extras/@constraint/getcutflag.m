function c = getcutflag(X)

c = [];
for i = 1:length(X.clauses)
    c = [c;X.clauses{i}.cut];
end