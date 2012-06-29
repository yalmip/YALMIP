function y = le(X,Y)
%LE (overloaded)

% Author Johan Löfberg
% $Id: le.m,v 1.3 2005-06-17 13:02:01 joloef Exp $

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
