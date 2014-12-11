function Constraint = dd(X)

if issymmetric(X)
    Constraint = [diag(X) >= sum(abs(X-diag(diag(X))),2)];
else
    error('dd requires a symmetric argument.');
end