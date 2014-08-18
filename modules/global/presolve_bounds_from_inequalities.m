function p = presolve_bounds_from_inequalities(p)
if p.K.l >0
    nnz_per_row = (p.F_struc | p.F_struc)*[0;ones(size(p.F_struc,2)-1,1)];
    valid_rows = find(nnz_per_row>1);
    valid_rows(valid_rows<=p.K.f)=[];
    valid_rows(valid_rows>p.K.f + p.K.l)=[];
    for j = valid_rows(:)'
        b = p.F_struc(j,1);
        a = p.F_struc(j,2:end);
        if nnz(p.F_struc(j,2:end))>1
            ap = a.*(a>0);
            am = a.*(a<0);
            for k = find(a)
                L = p.lb;
                U = p.ub;
                L(k) = 0;
                U(k) = 0;
                if a(k) > 0 & (p.ub(k)-p.lb(k)) > 1e-8
                    newlower = (-b - ap*U - am*L)/a(k);
                    p.lb(k) = max(p.lb(k),newlower);
                elseif a(k) < 0
                    newupper = (b + ap*U + am*L)/(-a(k));
                    p.ub(k) = min(p.ub(k),newupper);
                end
            end
        end
    end
end
