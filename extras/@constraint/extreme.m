function V = extreme(P)
% extreme Computes vertices (very rudimentary implementation)
%
% V = extreme(F)
%
% F : Constraint object defining a polytope
%
% See also CHEBYBALL, BOUNDINGBOX

V = extreme(lmi(P));