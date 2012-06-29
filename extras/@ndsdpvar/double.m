function X = double(X)
% DOUBLE (overloaded)

% Author Johan Löfberg
% $Id: double.m,v 1.2 2006-07-13 19:40:59 joloef Exp $

X = reshape(double(sdpvar(X)),X.dim);