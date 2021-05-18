function V = vertex(P)
% VERTEX Computes vertices (very rudimentary implementation)
%
% V = vertex(P)
%
% P : Constraint object defining a polytope
%
% See also CHEBYBALL, BOUNDINGBOX

V = vertex(lmi(P));