function [p,infeasible] = detectRedundantInfeasibleEXPRows(p)
infeasible = 0;
if any(p.K.e)
    newinEqualities = [];
    top = startofEXPCone(p.K);
    for j = 1:p.K.e
        % Looking for exp(x1/0)*0 <= x3
        m = p.K.e(j);
        M = p.F_struc(top:top+3-1,:);
        fixed = ~any(M(2,2:end));
        if fixed && M(2,1) == 0
            inequalityRows = [-M(1,:);M(3,:)]
            newinEqualities = [newinEqualities;inequalityRows];
            p.F_struc(top:top+3-1,:) = [];
        else
            top = top + 3;
        end
    end
    p = addInequality(p,newinEqualities);        
    p.K.e = p.K.e - size(newinEqualities,1)/2;
end