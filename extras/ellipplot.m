function [xe,ye]=ellipplot(P,gamma,ecolor,xc)
%DUAL2CELL Internal function for plotting ellipsoid

% Author Johan Löfberg
% $Id: ellipplot.m,v 1.3 2008-05-16 12:55:27 joloef Exp $

if nargin<4
   xc = zeros(length(P),1);
end

if nargin<3
   ecolor = [1 0 0];
end

if nargin<2
   gamma = 1;
end

P=P/gamma;
P=chol(P);

theta = linspace(-pi,pi,1000);
z = [cos(theta); sin(theta)];

x = inv(P)*z;
for n=1:length(x)
  x(:,n)=x(:,n)+xc;
end;

fill(x(1,:),x(2,:),ecolor)