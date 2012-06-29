function F = sdpvar(X)
% Internal class for constraint list

% Author Johan Löfberg
% $Id: sdpvar.m,v 1.1 2004-06-17 08:40:03 johanl Exp $


F = sdpvar(lmi(X));%X.Evaluated{1};
