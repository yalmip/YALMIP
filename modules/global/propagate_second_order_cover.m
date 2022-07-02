function [L,U] = propagate_second_order_cover(p_rawcuts,p)
L = -inf(length(p.c),1);
U = inf(length(p.c),1);
if p_rawcuts.K.l>0 && ~isinf(p.globalcardinality.up)
    for i = 1:p_rawcuts.K.l
        row = p_rawcuts.F_struc(i,:);
        [L_,U_] = presolve_glover_sherali_cover(row,p);
        if ~isempty(L_)
            L = max(L,L_);
            U = min(U,U_);
        end
    end
end