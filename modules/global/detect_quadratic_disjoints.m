function p = detect_quadratic_disjoints(p)

p.quadraticdisjoints = [];
if any(p.K.l)
    top = startofLPCone(p.K);
    for i = 1:p.K.l
        row = p.F_struc(top,:);
        if nnz(row)<=3
            a = row(1);            
            [~,var,val] = find(row(2:end));
            if length(var)==1
                if p.variabletype(var)==2
                    if val > 0 && a < 0
                        x = find(p.monomtable(var,:));                        
                        p.quadraticdisjoints = [p.quadraticdisjoints [x;sqrt(-a);-sqrt(-a)]];
                        % val*x^2 - a >= 0
                        % i.e >=sqrt(a) or <=-sqrt(a)
                    end
                end
            end
        end
        top = top + 1;
    end
end
