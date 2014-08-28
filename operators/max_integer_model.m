function F = max_integer_model(X,t);

[M,m] = derivebounds(X);

if all(M==m)
    F = [t == max(M)];
    return
end

n = length(X);
d = binvar(n,1);
F = (sum(d)==1);
F = F + (-(max(M)-min(m))*(1-d) <= t-X <= (max(M)-min(m))*(1-d));
kk = [];
ii = [];
for i = 1:n
    k = [1:1:i-1 i+1:1:n]';
    ii = [ii;repmat(i,n-1,1)];
    kk = [kk;k];
    Mm = M(k)-m(i);
end
xii = extsubsref(X,ii);
dii = extsubsref(d,ii);
xkk = extsubsref(X,kk);
F = F + (xkk <= xii+(M(kk)-m(ii)).*(1-dii));