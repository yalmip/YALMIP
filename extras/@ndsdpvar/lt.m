function y = lt(X,Y)
%LT (overloaded)

% Author Johan Löfberg
% $Id: lt.m,v 1.2 2006-07-13 19:40:59 joloef Exp $

if isa(X,'ndsdpvar')
    X = sdpvar(X);
elseif isa(X,'double')
    X = X(:);
end

if isa(Y,'ndsdpvar')
    Y = sdpvar(Y);
elseif isa(Y,'double')
    Y = Y(:);
end

try
    y = constraint(X,'<',Y);
catch
    error(lasterr)
end
