function y = ge(X,Y)
%GE (overloaded)

try
    y = constraint(X,'>=',Y);
catch
    error(lasterr)
end
