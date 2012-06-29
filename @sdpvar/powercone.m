function y = cone(z,x,y,alpha)
%CONE Defines a power cone x^alpha y ^(1-alpha) > |z|
%
% Input
%    z,y,x   : sclar SDPVAR objects.
%    alpha   : scalar double 0<=alpha<=1
%
% Example
%    F = set(cone(z,x,y,alpha))
%
% An alternative syntax with only one argument is also possible
%    F = set(cone(z))
% This command is equivalent to set(cone(z(2:end),z(1))
%

%
% See also  SET, CONE

% Author Johan Löfberg
% $Id: powercone.m,v 1.1 2008-06-27 13:07:46 joloef Exp $

if numel(z)>1
    error('x must be a scalar')
end
if numel(x)>1
    error('x must be a scalar')
end
if numel(y)>1
    error('y must be a scalar')
end

try
    y = [x;y;z;alpha];
    y.typeflag = 20;
    y = set(y);
catch
    rethrow(lasterror)
end