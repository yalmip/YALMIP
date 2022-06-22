function p = addInequality(p,row)
if ~anyCones(p.K)
    % Just append in the end
    p.F_struc = [p.F_struc;row];
else
    % Insert before conics
    if p.K.f==0 % a bit faster, append on top
        p.F_struc = [row;
                     p.F_struc];
    else
        p.F_struc = [p.F_struc(1:(p.K.f+p.K.l),:);
                     row;
                     p.F_struc(1+p.K.f+p.K.l:end,:)];
    end
end
p.K.l = p.K.l + size(row,1);