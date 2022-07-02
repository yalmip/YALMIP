function p = propagate_cardinality(p)

if p.feasible
    % Look in UBs and cardinality constraints
    for j = find(ismember(p.knapsack.type,[2 3 4 5]))
        xy = p.knapsack.variables{j};
        locked = p.lb(xy) == 1;
        if any(locked)
           % locked = find(p.lb(xy) == 1);
            m = nnz(locked);
            if m > p.knapsack.b{j}
                p.feasible = 0;
                return
            elseif m == p.knapsack.b{j}
                p.ub(p.knapsack.variables{j}) = 0;
                p.ub(xy(locked))=1;
                p.lb(xy(locked))=1;
            end
        end
    end
end