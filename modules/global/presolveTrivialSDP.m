function p = presolveTrivialSDP(p)

x0 = ones(length(p.c),1);
x0(find(p.lb == 0 & p.ub==0)) = 0;
if any(x0 == 0)
    newEqualities = [];
    if any(p.K.s)
        top = startofSDPCone(p.K);        
        for j = 1:length(p.K.s)
            X = p.F_struc(top:top + p.K.s(j)^2-1,:);
            X = X | X;
            X = reshape(X*[1;x0],p.K.s(j),p.K.s(j));
            % Look for zero diagonal. This means we can move all
            % nonzero elements to a zero equality
            e = find(diag(X)==0);
            if any(e)
                Z = spalloc(p.K.s(j),p.K.s(j),length(e)*2*p.K.s(j));
                for k = 1:length(e)
                    Z(:,e(k))=1;
                    Z(e(k),:)=1;
                end
                m1 = find(Z(:));      % To be removed
                m2 = find(triu(Z,1)); % To be moved
                equalityRows = p.F_struc(top + m2 - 1,:);
                p.F_struc(top + m1 - 1,:) = [];
                p.K.s(j) = p.K.s(j) - length(e);
                equalityRows = equalityRows(find(any(equalityRows,2)),:);
                newEqualities = [newEqualities;equalityRows];
            end
             top = top + p.K.s(j)^2;
        end
    end
    if ~isempty(newEqualities)
        keep = find(any(newEqualities,2));
        p.F_struc = [newEqualities(keep,:);p.F_struc];
        p.K.f = p.K.f + length(keep);
    end
end