function [Ax,Ay,b,K] = convexhullFromSampled(z,fz,xL,xU);

k1 = max((fz(2:end)-fz(1))./(z(2:end)-xL))+1e-12;
k2 = min((fz(2:end)-fz(1))./(z(2:end)-xL))-1e-12;
k3 = min((fz(1:end-1)-fz(end))./(z(1:end-1)-xU))+1e-12;
k4 = max((fz(1:end-1)-fz(end))./(z(1:end-1)-xU))-1e-12;
Ax = [-k1;k2;-k3;k4];
Ay = [1;-1;1;-1];
b =  [k1*(-z(1)) + fz(1);-(k2*(-z(1)) + fz(1));k3*(-z(end)) + fz(end);-(k4*(-z(end)) + fz(end))];
K.f = 0;
K.l = length(b);