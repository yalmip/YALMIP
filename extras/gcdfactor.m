function d = gcdfactor(A)
d = A(1);
if d~=round(d)
    d = 1;return
end
for n = 2:numel(A)
    if d == 1 
        return
    elseif A(n)~=round(A(n))
        d = 1;
        return
    else
        d = gcd(d,A(n));
    end
end