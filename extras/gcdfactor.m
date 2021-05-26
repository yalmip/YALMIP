function d = gcdfactor(A)
d = A(1);
for n = 2:numel(A)
    if d == 1
        return
    else
        d = gcd(d,A(n));
    end
end