function y = cone(Axplusb,cxplusd)
%CONE Defines a second order cone constraint ||z||<x
%
% Input
%    z       : SDPVAR object.
%    h       : SDPVAR object
%
% Example
%
% Standard SOCP constraint ||z||<x where z is a vector and x is a scalar
%    F = set(cone(z,x)) 
%
% Alternative syntax with only one argument is also possible
%    F = set(cone(z))
% This command is equivalent to set(cone(z(2:end),z(1))
%
% To quickly define several cones, the argument can be a matrix, and the
% command is then short-hand for 
% for i = 1:size(z,2);F = [F,cone(z(:,i))];end 
% The code will however run much fast than the manual version
%
% See also  SET, RCONE, @SDPVAR/NORM

% Author Johan Löfberg
% $Id: cone.m,v 1.7 2008-04-24 11:15:13 joloef Exp $

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
    y = set(y);
catch
    error(lasterr)
end