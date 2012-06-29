function v = depends(X)
% depends (overloaded)

% Author Johan Löfberg
% $Id: depends.m,v 1.2 2006-07-13 19:40:59 joloef Exp $

v = depends(sdpvar(X));