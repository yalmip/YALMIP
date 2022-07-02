function p = propagate_atleast(p)

if p.feasible
    for j = 1:length(p.atleast.groups)
        xy = p.atleast.groups{j};
        locked = p.ub(xy) == 0;
        if any(locked)
           % locked = find(p.lb(xy) == 1);
            m = nnz(locked);
            free = length(xy)-m;
            if free < p.atleast.bounds(j)
                p.feasible = 0;
                return
            elseif free == p.atleast.bounds(j)
                p.lb(p.atleast.groups{j}) = 1;
                p.ub(xy(locked))=0;
                p.lb(xy(locked))=0;
            end
        end
    end
end