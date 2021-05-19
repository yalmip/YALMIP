function [p,infeasible] = detectRedundantInfeasibleSOCPRows(p)
infeasible = 0;
if any(p.K.q)
    newEqualities = [];
    top = startofSOCPCone(p.K);
    for j = 1:length(p.K.q)
        % Looking for norm(x) <= k and k=0 or k<0
        m = p.K.q(j);
        M = p.F_struc(top:top+m-1,:);
        fixed = ~any(M(1,2:end));
        if fixed
            if M(1) <= -p.options.bnb.feastol
                infeasible = 1;
                return
            elseif M(1) <= 0
                % Not numerically infeasible,
                % but at least says x = 0
                equalityRows = M(2:m,:);
                newEqualities = [newEqualities;equalityRows];
                p.F_struc(top:top+m-1,:) = [];
                p.K.q(j) = 0;
            end
            top = top + p.K.q(j);
        end
    end
    p = addEquality(p,newEqualities);        
    p.K.q(p.K.q == 0) = [];
end