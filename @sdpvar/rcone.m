function y = rcone(z,x,y)
%RCONE Defines a rotated second order cone constraint ||z||^2<2xy, x+y>0
%
% Input
%    z       : column vector SDPVAR object.
%    h       : scalar double or SDPVAR object
%
% Example
%
%    F = (rcone(z,x,y)) 
%
% See also CONE


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
  y = lmi(y);
catch
  error(lasterr)
end