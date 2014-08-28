function c = getcutflag(X)

X = flatten(X);
c = [];
for i = 1:length(X.clauses)
    c = [c;X.clauses{i}.cut];
end