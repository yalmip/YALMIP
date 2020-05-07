function y = pcone(z,x,alpha)
%PCONE Defines a power cone x(1)^alpha x(2)^(1-alpha) > |z|
%
% Input
%    z,x     : SDPVAR objects.
%    alpha   : scalar double 0<=alpha<=1
%
% Example
%    F = pcone(z,x,alpha)
%
% An alternative syntax with only one argument is also possible
%    F = pcone(Z) # where Z = [x;z;alpha]
%
%
% See also @sdpvar/CONE, @sdpvar/EXPCONE

try
    y = [reshape(x,[],1);reshape(z,[],1);reshape(alpha,[],1)];
    y.typeflag = 20;
    y = lmi(y);
catch
    rethrow(lasterror)
end