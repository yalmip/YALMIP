function [p,infeasible] = detectRedundantInfeasibleSDPRows(p)
infeasible = 0;
if any(p.K.s)
    top = startofSDPCone(p.K);
    newEqualities = [];
    for j = 1:length(p.K.s)
        X = p.F_struc(top:top + p.K.s(j)^2-1,:);
        X = reshape(any(X,2),p.K.s(j),p.K.s(j));
        e = find(~any(X,2));
        if any(e)
            % Not a single element is used, so simply and reduce SDP
            Z = spalloc(p.K.s(j),p.K.s(j),length(e)*2*p.K.s(j));
            for k = 1:length(e)
                Z(:,e(k))=1;
                Z(e(k),:)=1;
            end
            m = find(Z(:));
            p.F_struc(top + m - 1,:)=[];
            p.K.s(j) = p.K.s(j) - length(e);
        else
            % Look for zero diagonal. This means we can move all
            % nonzero elements to a zero equality
            e = find(diag(X)==0);
            if any(e)
                Z = spalloc(p.K.s(j),p.K.s(j),length(e)*2*p.K.s(j));
                for k = 1:length(e)
                    Z(:,e(k)) = 1;
                    Z(e(k),:) = 1;
                end
                m1 = find(Z(:)); % To be removed
                m2 = find(triu(Z,1)); % To be moved
                equalityRows = p.F_struc(top + m2 - 1,:);
                p.F_struc(top + m1 - 1,:) = [];
                p.K.s(j) = p.K.s(j) - length(e);
                equalityRows = equalityRows(find(any(equalityRows,2)),:);
                newEqualities = [newEqualities;equalityRows];
            end
        end
        top = top + p.K.s(j)^2;
    end
    if ~isempty(newEqualities)
        p = addEquality(p,newEqualities);
    end
    p.K.s(p.K.s == 0) = [];
end
