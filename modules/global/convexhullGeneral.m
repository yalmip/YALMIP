function [Ax,Ay,b,K] = convexhullGeneral(xL,xU,f)

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

[Ax,Ay,b] = convexhullFromSampled(z,fz,xL,xU);