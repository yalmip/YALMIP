function t = alldifferent(X)
% alldifferent (overloaded)

% Author Johan Löfberg
% $Id: alldifferent.m,v 1.1 2006-05-17 15:15:25 joloef Exp $

t = alldifferent(sdpvar(X));