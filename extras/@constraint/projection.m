function  F = projection(F,x,varargin)
% projection  Projects polytopic set object (Requires the Multi-parametric Toolbox).
%
% Fproj     = projection(F,x)
%
% F      : Polytopic set object
% x      : Variables to project on
% method : See HELP PROJECTION

% Author Johan Löfberg
% $Id: projection.m,v 1.1 2010-03-01 15:22:58 joloef Exp $

F = projection(set(F),x,varargin{:});
