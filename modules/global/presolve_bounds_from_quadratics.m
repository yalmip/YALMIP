function p = presolve_bounds_from_quadratics(p)
if p.K.l >0
    if any(p.variabletype==2)
        for i = 1:p.K.l
            U = p.F_struc(p.K.f+i,1);
            if U>=0
                [aux,col,val] = find(p.F_struc(p.K.f+i,2:end));
                col = col(:);
                if all(p.variabletype(col)==2) & all(val<=0)
                    p.ub(col) = min([p.ub(col) U./-val(:)],[],2);
                    p.lb(col) = max([p.lb(col) 0*val(:)],[],2);                 
                end
            end
        end
    end
end