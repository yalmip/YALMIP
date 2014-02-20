function y = le(X,Y)
%LE (overloaded)

try
    y = constraint(X,'<=',Y);
catch
    error(lasterr)
end
