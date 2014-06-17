function [Ax,Ay,b,K] = convexhullConvex(varargin)
% Two upper bounds from tangents
% y > f(xL) + (x-xL)*df(xL)
% y > f(xU) + (x-xL)*df(xU)
% Upper bound from conneting extreme points
% y < f(xU)(x-xL)/(xU-xL) +  f(xL)(xU-x)/(xU-xL)
% can be wrtitten as
% Ax*x + Ay*y < b

if rem(nargin,3)
    error('The convex hull generator assumes n triplets (x,f,df).')
end
m = nargin/3;
x = [varargin{(1:m)}]';
f = [varargin{(1:m)+m}]';
df = [varargin{(1:m)+2*m}]';

if all(diff(x)>0)
    Ay = [-ones(m,1);1];
    b  = [-f + x.*df; -f(end)*x(1)/(x(end)-x(1)) +  f(1)*x(end)/(x(end)-x(1))];
    Ax  = [df;-f(end)/(x(end)-x(1)) + f(1)/(x(end)-x(1))];
else
    Ax = [];
    Ay = [];
    b = [];
end
K.f = 0;
K.l = length(b);