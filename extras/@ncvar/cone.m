function y = cone(Axplusb,cxplusd)
%CONE Defines a second order cone constraint ||z||<x
%
% Input
%    z       : column vector SDPVAR object.
%    h       : scalar double or SDPVAR object
%
% Example
%    F = cone(z,x))
%
% An alternative syntax with only one argument is also possible
%    F = cone(z)
% This command is equivalent to cone(z(2:end),z(1)
%
% See also RCONE, @SDPVAR/NORM

[n,m] = size(Axplusb);
if min(n,m)>1
    error('z must be a  vector')
end

if nargin == 2
    if prod(size(cxplusd))>1
        error('x must be a scalar')
    end
else
end

if n<m
    Axplusb = Axplusb';
end

try
    if nargin == 2
        y = [cxplusd;Axplusb];
    else
        y = [Axplusb];
    end
    if is(y,'linear')
        y.typeflag = 4;
    else
        y = (Axplusb'*Axplusb < cxplusd*cxplusd);
    end
catch
    error(lasterr)
end