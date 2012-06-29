function v = isreal(X)
% isreal (overloaded)

% Author Johan Löfberg
% $Id: isreal.m,v 1.2 2006-07-25 12:57:08 joloef Exp $

v = isreal(sdpvar(X));