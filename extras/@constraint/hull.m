function [Fhull,t,y] = hull(varargin)
% HULL  Construct a model of the convex hull
%
% H = hull(F1,F2,...)
%
% OUTPUT
%   H   : SET object describing the convex hull of the input constraints
%
% INPUT
%   Fi  : SET objects with constraints
%
% Note that the convex representation of the convex hull requires a lifting
% (introduction of new variables). Hence, if you have many set of
% constraints, your problem rapidly grows large.

% $Id: hull.m,v 1.2 2008-02-14 14:53:36 joloef Exp $   

for i = 1:nargin  
    if isa(varargin{i},'constraint')
        varargin{i} = set(varargin{i});
    end
end
% Call the hull operator from the set
[Fhull,t,y] = hull(varargin{:});