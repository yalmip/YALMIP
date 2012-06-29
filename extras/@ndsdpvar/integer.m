function x = integer(x)
% integer (overloaded)

% Author Johan Löfberg
% $Id: integer.m,v 1.1 2006-05-17 12:03:46 joloef Exp $

x = sdpvar(x);
x = binary(x);