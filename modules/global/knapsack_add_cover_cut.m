function cuts = knapsack_add_cover_cut(p,x,alg,upper)

if ~isinf(upper) && p.LinearBinaryPositiveCost
    % We can add a global knapsack from c'*x <= upper
    objective_upper = upper - p.f;
    if objective_upper >= 0
        relevant_to_use = find(p.c & p.c <= objective_upper);
        if ~isempty(relevant_to_use)
            p.knapsack.a{end+1} = p.c(relevant_to_use);
            p.knapsack.b{end+1} = objective_upper;
            p.knapsack.variables{end+1} = relevant_to_use;
            p.knapsack.type(end+1) = 0;
        end
    end
end

cuts = [];
for i = find(p.knapsack.type == 0)
    a = p.knapsack.a{i};
    b = p.knapsack.b{i};
    v = p.knapsack.variables{i};
    x_ = x(v);
    
    cut_ = knapsack_create_cover_cut(a,b,x_,alg);
    if ~isempty(cut_)
        cut = spalloc(1,length(p.c)+1,0);
        cut(1) = cut_(1);
        cut(1 + v) = cut_(2:end);
        cuts = [cuts;cut];
    end
end

