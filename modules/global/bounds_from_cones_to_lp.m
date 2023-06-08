function p = bounds_from_cones_to_lp(p)

if any(p.K.e)
    % x3 >= x2*exp(x1/x2) so x2>=0
    top = startofEXPCone(p.K);
    newInequalities = [];
    for i = 1:p.K.e
        newInequalities = [newInequalities;p.F_struc(top+1,:)];
        top = top + 3;
    end
    p = addInequality(p,newInequalities);
end

if any(p.K.p)
    % x1^a x2^(1-a)>=norm(x3) so x1,x2>=0
    top = startofPOWCone(p.K);
    newInequalities = [];
    for i = 1:length(p.K.p)
        newInequalities = [newInequalities;p.F_struc(top:top+1,:)];
        top = top + p.K.p(i);
    end
    p = addInequality(p,newInequalities);
end

if any(p.K.s)
    % Diagonals are non-negative
    top = startofSDPCone(p.K);
    newInequalities = [];
    for i = 1:length(p.K.s)
        n = p.K.s(i);
        index = top + find(speye(n))-1;        
        row = p.F_struc(index,:);
        for k = 1:n
            if nnz(row(k,2:end))> 1
                % newInequalities = [newInequalities;row];        
            elseif nnz(row(k,2:end)) == 1
                % This is just a bound
                idx = find(row(k,2:end));
                if row(k,idx+1)>0
                    p.lb(idx) = max(p.lb(idx),-row(k,1)/row(k,idx+1));
                else
                    p.ub(idx) = min(p.ub(idx),row(k,1)/(-row(k,idx+1)));
                end
            end
        end
        top = top + n^2;
    end
    newInequalities = unique(newInequalities,'rows');
    relevant = any(newInequalities(:,2:end),2);
    newInequalities = newInequalities(relevant,:);
    p = addInequality(p,newInequalities);
end