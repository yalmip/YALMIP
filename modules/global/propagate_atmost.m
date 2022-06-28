function p = propagate_atmost(p)

if p.feasible
    for j = 1:length(p.atmost.groups)
        xy = p.atmost.groups{j};
        locked = p.lb(xy) == 1;
        if any(locked)
           % locked = find(p.lb(xy) == 1);
            m = nnz(locked);
            if m > p.atmost.bounds(j)
                p.feasible = 0;
                return
            elseif m == p.atmost.bounds(j)
                p.ub(p.atmost.groups{j}) = 0;
                p.ub(xy(locked))=1;
                p.lb(xy(locked))=1;
            end
        end
    end
end