function y = lt(X,Y)
%LT (overloaded)

% Author Johan Löfberg
% $Id: lt.m,v 1.1 2006-08-10 18:00:21 joloef Exp $

if isa(X,'blkvar')
    X = sdpvar(X);
end

if isa(Y,'blkvar')
    Y = sdpvar(Y);
end

try
    y = constraint(X,'<',Y);
catch
    error(lasterr)
end
