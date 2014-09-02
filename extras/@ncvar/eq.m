function y = eq(X,Y)
%EQ (overloaded)

if isa(X,'blkvar')
    X = sdpvar(X);
end

if isa(Y,'blkvar')
    Y = sdpvar(Y);
end

try
    y = constraint(X,'==',Y);
catch
    error(lasterr)
end
