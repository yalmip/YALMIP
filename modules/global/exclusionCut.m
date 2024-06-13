function [b,a] = exclusionCut(x,thesign)
% Returns cut [b a]*[1;x] which is infeasible 
% for given binary
% thesign = -1 models same for negated binary

zv = find(abs(x) <= 1e-5);
nz = find(abs(x) > 1e-5);   
a = spalloc(1,length(x),length(zv));
b = length(x)-length(zv)-1;
a = sparse(1,zv,1,1,length(x))-sparse(1,nz,1,1,length(x));
if thesign == -1
    a = -a;
end