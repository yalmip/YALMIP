function y = rcone(z,x,y)
%RCONE Defines a rotated second order cone constraint ||z||^2<2xy, x+y>0
%
% Input
%    z       : column vector SDPVAR object.
%    h       : scalar double or SDPVAR object
%
% Example
%
%    F = set(rcone(z,x,y)) 
%
% See also  SET, CONE


% Author Johan Löfberg
% $Id: rcone.m,v 1.1 2006-08-10 18:00:22 joloef Exp $

if size(z,2)>1
  error('z must be a column vector')
end

if prod(size(x))>1
  error('x must be a scalar')
end
if prod(size(y))>1
  error('y must be a scalar')
end

try
  y = [x;y;z];
  y.typeflag = 5;
catch
  error(lasterr)
end