function x = value(x);
%VALUE (Overloaded)

x = value(x.Evaluated{1});

if isessentiallysymmetric(x)
    x = min(eig(x)) >= 0;
else
    x = min(x(:)) >= 0;
end
