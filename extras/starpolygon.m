function F = starpolygon(xi,yi,z,c,t,p)
%PWF Defines a star-shaped polygon set
%
% F = STARPOLYGON(xi,yi,z,c,t,p)
%
% xi,yi : Coordinates (DOUBLE)
% z     : Variable (2x1 SDPVAR)
% c     : Translation (2x1 SDPVAR/DOUBLE, optional)
% t     : Scale (Scalar DOUBLE/SDPVAR, optional)
% p     : Vantage point (2x1 DOUBLE, optional)
%
% Default is c = [0;0], t = 1, p = [0;0], i.e. it is assumed that the data
% represents a set star-convex around the origin. 

if nargin < 4 || isempty(c)
    c = [0;0];
end
if nargin < 5 || isempty(t)
    t = 1;
end
if nargin < 6 || isempty(p)
    p = [0;0];
elseif isa(p,'double')
else
    xc = mean(xi);xi = xi-xc;
    yc = mean(yi);yi = yi-yc;
end
lambda = sdpvar(length(xi),1);
F = [sos2(lambda),lambda>=0
    z == c + p + (lambda'*[xi(:)-p(1) yi(:)-p(2)])'];
F = [F, sum(lambda)<=t];
