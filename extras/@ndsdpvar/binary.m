function x = binary(x)
% binary (overloaded)

% Author Johan Löfberg
% $Id: binary.m,v 1.2 2006-07-13 19:40:58 joloef Exp $

x = binary(sdpvar(x));