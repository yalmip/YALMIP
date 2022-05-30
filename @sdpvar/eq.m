function y = eq(X,Y)
%EQ (overloaded)

if isa(X,'blkvar')
    X = sdpvar(X);
end

if isa(Y,'blkvar')
    Y = sdpvar(Y);
end

try
    if ishermitian(X) && ishermitian(Y)
        y = constraint(triu(X),'==',triu(Y));
    else
        y = constraint(X,'==',Y);
    end
catch
    error(lasterr)
end
