function [Ax,Ay,b,K] = convexhullConvex(xL,xU,f)

z = linspace(xL,xU,100);
fz = f(z);
[minval,minpos] = min(fz);
[maxval,maxpos] = max(fz);
xtestmin = linspace(z(max([1 minpos-5])),z(min([100 minpos+5])),100);
xtestmax = linspace(z(max([1 maxpos-5])),z(min([100 maxpos+5])),100);

fz1 = f(xtestmin);
fz2 = f(xtestmax);
z = [z(:);xtestmin(:);xtestmax(:)];
fz = [fz(:);fz1(:);fz2(:)];
[z,sorter] = sort(z);
fz = fz(sorter);
[z,ii,jj]=unique(z);
fz = fz(ii);

k1 = max((fz(2:end)-fz(1))./(z(2:end)-xL))+1e-12;
k2 = min((fz(2:end)-fz(1))./(z(2:end)-xL))-1e-12;
k3 = min((fz(1:end-1)-fz(end))./(z(1:end-1)-xU))+1e-12;
k4 = max((fz(1:end-1)-fz(end))./(z(1:end-1)-xU))-1e-12;
Ax = [-k1;k2;-k3;k4];
Ay = [1;-1;1;-1];
b =  [k1*(-z(1)) + fz(1);-(k2*(-z(1)) + fz(1));k3*(-z(end)) + fz(end);-(k4*(-z(end)) + fz(end))];
