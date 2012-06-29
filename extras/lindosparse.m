function [Nbegcol,Nlencol,Nrowndx] = lindosparse(cJacobian);

Nbegcol = [];
Nrowndx = [];
Nlencol = [];
top = 0;
for i = 1:size(cJacobian,2)
    [ii,jj,kk] = find(cJacobian(:,i));
    if isempty(ii)
        Nbegcol = [Nbegcol top];
        Nlencol = [Nlencol 0];
    else
        Nbegcol = [Nbegcol top];
        Nrowndx = [Nrowndx ii(:)'-1];
        Nlencol = [Nlencol length(ii)];
        top = top + length(ii);
    end
end
if  isempty(Nrowndx)
    Nrowndx = [];
end
Nbegcol = [Nbegcol sum(Nlencol)];