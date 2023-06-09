function V = vertex(varargin)
% VERTEX Computes vertices (very rudimentary implementation)
%
% V = VERTEX(P)
%
% P : Constraint object defining a polytope
%
% With a second argument x the projection of the vertices
% to the x-space is returned.
%
% See also CHEBYBALL, BOUNDINGBOX, POLYTOPE

varargin{1} = lmi(varargin{1});
V = vertex(varargin{:});