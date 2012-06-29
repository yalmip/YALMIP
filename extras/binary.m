function x = binary(x)
%BINARY Overloaded
%
% Author Johan Löfberg
% $Id: binary.m,v 1.1 2005-12-01 14:09:48 joloef Exp $

if isempty(x)
    x = [];
else
    error('BINARY can only be applied to SDPVAR objects or empty doubles');
end