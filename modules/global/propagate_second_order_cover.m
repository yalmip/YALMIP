function [L,U] = propagate_second_order_cover(cuts,p)
L = -inf(length(p.c),1);
U = inf(length(p.c),1);
if cuts.K.l>0 && ~isinf(p.binarycardinality.up)
    for i = 1:cuts.K.l
        row = cuts.F_struc(cuts.K.f+i,:);
        [L_,U_] = presolve_glover_sherali_cover(row,p);
        if ~isempty(L_)
            if any(L_>L)
                L = max(L,L_);
            end
            if any(U_<U)
                U = min(U,U_);
            end
        end
    end
end