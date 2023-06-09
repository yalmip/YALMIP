function p = presolve_remove_equalitydominateinequality(p)
% Look for b == a'*x, c>=a'*x where c>b
if p.K.f > 0 && p.K.l > 0 && p.feasible
    remove = [];
    hash = p.F_struc(1:p.K.f+p.K.l,:)*[0;randn(length(p.c),1)];
    [~,idx,mappedto] = unique(hash,'stable');
    if length(idx) < p.K.f+p.K.l
        for i = 1:length(idx)
            if idx(i) <= p.K.f
                j = find(mappedto == idx(i));
                if length(j)>1
                    j = j(j > p.K.f);
                    redundant = j(find(p.F_struc(j,1) > p.F_struc(idx(i),1)));
                    remove = [remove redundant];
                end
            end
        end
        if ~isempty(remove)
            remove = unique(remove);
            p.F_struc(remove,:)=[];
            p.K.l = p.K.l - length(remove);
        end
    end
end