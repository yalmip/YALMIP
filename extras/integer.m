function x = integer(x)
%INTEGER Overloaded
%
% Author Johan Löfberg
% $Id: integer.m,v 1.1 2005-12-01 14:09:48 joloef Exp $

if isempty(x)
    x = [];
else
    error('INTEGER can only be applied to SDPVAR objects or empty doubles');
end