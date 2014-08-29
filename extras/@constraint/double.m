function x = double(x);
%DOUBLE (Overloaded)

x = double(x.Evaluated{1});

if issymmetric(x)
    x = min(eig(x)) >= 0;
else
    x = min(x(:)) >= 0;
end
