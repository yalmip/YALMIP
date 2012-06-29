function X = power(X,d)
% UMINUS (overloaded)

% Author Johan Löfberg
% $Id: power.m,v 1.1 2008-04-24 08:27:58 joloef Exp $

s = X.dim;
X = sdpvar(X);
d = d(:);
X = power(X,d);
X = reshape(X,s);