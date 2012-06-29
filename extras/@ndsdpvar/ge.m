function y = ge(X,Y)
%GE (overloaded)

% Author Johan Löfberg
% $Id: ge.m,v 1.3 2006-07-25 12:57:08 joloef Exp $

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
    y = constraint(X,'>=',Y);
catch
    error(lasterr)
end
