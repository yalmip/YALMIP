function [p,infeasible] = detectRedundantInfeasiblePOWRows(p)
infeasible = 0;
if any(p.K.p)
    newEqualities = [];
    top = startofPOWCone(p.K);
    for j = 1:length(p.K.p)
        % Looking for 0 >= norm(x3)
        m = p.K.p(j);
        M = p.F_struc(top:top+1,:);
        fixed = ~any(M,2);
        if any(fixed)           
            newEqualities = [newEqualities;
                             p.F_struc(top+2:top+m-2,:)];
            p.F_struc(top:top+m-1,:) = [];
            p.K.p(j) = 0;
        else
            top = top + m;
        end
    end
    p = addEquality(p,newEqualities);        
    p.K.p(p.K.p == 0) = [];
end