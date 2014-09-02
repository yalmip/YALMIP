function y = lt(X,Y)
%LE (overloaded)

if isa(X,'blkvar')
    X = sdpvar(X);
end

if isa(Y,'blkvar')
    Y = sdpvar(Y);
end

try
    y = constraint(X,'<=',Y);
catch
    error(lasterr)
end
