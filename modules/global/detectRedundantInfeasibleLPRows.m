function [p,infeasible] = detectRedundantInfeasibleLPRows(p)
infeasible = 0;
if any(p.K.f) || any(p.K.l)
    zero_row_eq = find(~any(p.F_struc(1:p.K.f,2:end),2));
    zero_row_lp = find(~any(p.F_struc(1+p.K.f:p.K.f+p.K.l,2:end),2));
    if ~isempty(zero_row_lp)
        lhs = p.F_struc(p.K.f + zero_row_lp,1);
        zero_row_pos = find(lhs >= 0);
        remove_these = zero_row_lp(zero_row_pos);
        p.F_struc(p.K.f + remove_these,:) = [];
        p.K.l = p.K.l - length(remove_these);
        zero_row_neg = find(lhs < -p.options.bnb.feastol);
        if ~isempty(zero_row_neg)
            infeasible = 1;
            return
        end
    end
    if ~isempty(zero_row_eq)
        lhs = p.F_struc(zero_row_eq,1);
        zero_row_ok = find(abs(lhs)<p.options.bnb.feastol);
        remove_these = zero_row_lp(zero_row_ok);
        p.F_struc(remove_these,:) = [];
        p.K.f = p.K.f - length(remove_these);
        zero_row_bad = find(abs(lhs)>p.options.bnb.feastol);
        if ~isempty(zero_row_bad)
            infeasible = 1;
            return
        end
    end
end