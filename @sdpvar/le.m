function y = le(X,Y)
%LE (overloaded)

% Author Johan Löfberg

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
