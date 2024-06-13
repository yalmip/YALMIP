function y = cone(Axplusb,cxplusd)
%CONE Defines a second order cone constraint norm(z,2)<=x
%
% Input
%    z       : Linear SDPVAR object.
%    x       : Linear scalar SDPVAR object
%
% Example
%
% Standard SOCP constraint normally conveniently written as norm(z,2)<=x
% where z is a vector and x is a scalar 
%    F = cone(z,x)
%
% Alternative syntax with only one argument is also possible
%    F = cone(z)
% This is equivalent to cone(z(2:end),z(1))
%
% To quickly define several cones, the argument can be a matrix, and the
% command is then short-hand for 
% for i = 1:size(z,2);F = [F,cone(z(:,i))];end 
% The code will however run much faster than the non-vectorized version
%
% See also  @SDPVAR/EXPCONE, @SDPVAR/PCONE, @SDPVAR/NORM

[n,m] = size(Axplusb);
if min(n,m)>1 & nargin>1
    error('z must be a  vector')
end

if nargin == 2
    if prod(size(cxplusd))>1
        error('x must be a scalar')
    end
else
end

if n<m & nargin>1
    Axplusb = Axplusb';
end

if nargin > 1 & ~is(Axplusb,'real')
    y = cone([real(Axplusb);imag(Axplusb)],cxplusd);
    return
end

try
    if nargin == 2
        y = [cxplusd;Axplusb];
    else
        y = [Axplusb];
    end
    if size(Axplusb,2)>1
        y.typeflag = 54; % Stacked cones
    else
        y.typeflag = 4;
    end
    y = lmi(y);
catch
    error(lasterr)
end