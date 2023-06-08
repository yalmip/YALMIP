function p = propagate_setcover(p)

if p.feasible
    for j = 1:length(p.knapsack.b)
        if ismember(p.knapsack.type(j),[9 10])
            xy = p.knapsack.variables{j};
            locked = p.ub(xy) == 0;
            if any(locked)
                m = nnz(locked);
                free = length(xy)-m;
                if free < p.knapsack.b{j}
                    p.feasible = 0;
                    return
                elseif free == p.knapsack.b{j}
                    p.lb(p.knapsack.variables{j}) = 1;
                    p.ub(xy(locked))=0;
                    p.lb(xy(locked))=0;
                end
            end
        end
    end
end