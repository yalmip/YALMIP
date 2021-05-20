function y = powcone(x,y,z,alpha)
%POWCONE Defines a power cone x^alpha y^(1-alpha) > norm(z)
%
% Input
%    x,y,z     : SDPVAR objects or doubles
%    alpha     : scalar double 0<=alpha<=1
%
% Example
%    F = pcone(x,y,z,alpha)
%
% An alternative syntax with only one argument is also possible
%    F = pcone(Z) # where Z = [x;y;z;alpha]
%
% A vectorized version is also possible
%    F = pcone(Z) # where Z(:,i) = [x;y;z;alpha]
%
% See also @sdpvar/CONE, @sdpvar/EXPCONE

if nargin == 1
    [n,m] = size(x);
    if n < 4
        error('x must be a vector or matrix of height 4 or more ')
    end
    y = x;
    if m > 1
        y.typeflag = 58;
    else
        y.typeflag = 20;
    end
    y = lmi(y);
else
    y = powcone([x;y;reshape(z,[],1);alpha]);
end
