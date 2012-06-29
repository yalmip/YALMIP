function y = eq(X,Y)
%EQ (overloaded)

% Author Johan Löfberg
% $Id: eq.m,v 1.5 2006-05-17 13:22:51 joloef Exp $

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
