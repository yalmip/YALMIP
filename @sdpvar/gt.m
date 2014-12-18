function y = gt(X,Y)
%GT (overloaded)

if isa(X,'blkvar')
    X = sdpvar(X);
end

if isa(Y,'blkvar')
    Y = sdpvar(Y);
end

warning('YALMIP:strict','Strict inequalities are not supported. A non-strict has been added instead'')');

try
    y = constraint(X,'>=',Y);
catch
    error(lasterr)
end
