function [b,a] = exclusionCut(x,thesign)

zv = find(abs(x) <= 1e-5);
nz = find(abs(x) > 1e-5);   
a = spalloc(1,length(x),length(zv));
b = length(x)-length(zv)-1;
if thesign == -1
    a(nz) = 1;
    a(zv) = -1;
elseif thesign == 1        
    a(zv) = 1;
    a(nz) = -1;
end