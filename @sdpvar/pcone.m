function y = pcone(z,x,y,alpha)
%PCONE Defines a power cone x^alpha y ^(1-alpha) > |z|
%
% Input
%    z,y,x   : sclar SDPVAR objects.
%    alpha   : scalar double 0<=alpha<=1
%
% Example
%    F = set(pcone(z,x,y,alpha))
%
% An alternative syntax with only one argument is also possible
%    F = set(pcone(z))
%

%
% See also  SET, CONE

% Author Johan Löfberg
% $Id: pcone.m,v 1.2 2008-04-24 20:19:57 joloef Exp $

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