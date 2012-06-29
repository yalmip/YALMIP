function F = times(X,Y)
% times (overloaded)

% Author Johan Löfberg
% $Id: times.m,v 1.2 2006-07-28 06:27:01 joloef Exp $

if isa(X,'ndsdpvar')
    dim = X.dim;
    X = sdpvar(X);
else
    X = X(:);
end
if isa(Y,'ndsdpvar')
    dim = Y.dim;
    Y = sdpvar(Y);
else
    Y = Y(:);
end

F = X.*Y;
F = reshape(F,dim);