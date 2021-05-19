function p = addInequality(p,row)
if ~anyCones(p.K)
    % Just append in the end
    p.F_struc = [p.F_struc;row];
else
    % Insert before conics
    p.F_struc = [p.F_struc(1:(p.K.f+p.K.l),:);
                 row;
                 p.F_struc(1+p.K.f+p.K.l:end,:)];
end
p.K.l = p.K.l + size(row,1);