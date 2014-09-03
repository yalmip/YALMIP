function y = pcone(z,x,y,alpha)
%PCONE Defines a power cone x^alpha y ^(1-alpha) > |z|
%
% Input
%    z,y,x   : sclar SDPVAR objects.
%    alpha   : scalar double 0<=alpha<=1
%
% Example
%    F = pcone(z,x,y,alpha)
%
% An alternative syntax with only one argument is also possible
%    F = pcone(z)
%
%
% See also CONE

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
    y = lmi(y);
catch
    rethrow(lasterror)
end