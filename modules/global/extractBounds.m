function p = extractBounds(p)
if ~isempty(p.F_struc)
    [lb,ub,used_rows_eq,used_rows_lp] = find_lp_bounds(p.F_struc,p.K);
    if ~isempty(used_rows_lp)
        used_rows_lp = used_rows_lp(~any((p.F_struc(p.K.f + used_rows_lp,1+p.nonlinear)),2));
        if ~isempty(used_rows_lp)
            lower_defined = find(~isinf(lb));
            if ~isempty(lower_defined)
                p.lb(lower_defined) = max(p.lb(lower_defined),lb(lower_defined));
            end
            upper_defined = find(~isinf(ub));
            if ~isempty(upper_defined)
                p.ub(upper_defined) = min(p.ub(upper_defined),ub(upper_defined));
            end
            p.F_struc(p.K.f + used_rows_lp,:)=[];
            p.K.l = p.K.l - length(used_rows_lp);
        end
    end    
    if ~isempty(used_rows_eq)
        used_rows_eq = used_rows_eq(~any(full(p.F_struc(used_rows_eq,1+p.nonlinear)),2));
        if ~isempty(used_rows_eq)
            lower_defined = find(~isinf(lb));
            if ~isempty(lower_defined)
                p.lb(lower_defined) = max(p.lb(lower_defined),lb(lower_defined));
            end
            upper_defined = find(~isinf(ub));
            if ~isempty(upper_defined)
                p.ub(upper_defined) = min(p.ub(upper_defined),ub(upper_defined));
            end
            p.F_struc(used_rows_eq,:)=[];
            p.K.f = p.K.f - length(used_rows_eq);
        end
    end
end

