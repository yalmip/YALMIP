function y = powercone(z,x,y,alpha)
%POWERCONE Defines a power cone x^alpha y ^(1-alpha) > |z|
%
% Input
%    z,y,x   : sclar SDPVAR objects.
%    alpha   : scalar double 0<=alpha<=1
%
% Example
%    F = powercone(z,x,y,alpha)

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