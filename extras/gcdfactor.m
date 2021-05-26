function x = gcdfactor(a,b)
if nargin == 1
    if length(a)==1
        x = a;
        return
    else
        x = gcdfactor(a(1),a(2:end));
        return
    end
end
if length(b)==1
    x = gcd(a,b);
elseif all(a==b)
    x = a;
else
    x = gcd(a,b(1));
    if x == 1
        return
    else
        x = gcdfactor(x,b(2:end));
    end
end