function cuts = knapsack_add_cover_cut(p,x,alg,upper)

if ~isinf(upper) && p.LinearBinaryPositiveCost
    % We can add a global knapsack from c'*x <= upper
    relevant_to_use = find(p.c & p.c <= upper);
    p.knapsack.a{end+1} = p.c(relevant_to_use);
    p.knapsack.b{end+1} = upper;
    p.knapsack.variables{end+1} = relevant_to_use;
    p.knapsack.type(end+1) = 0;    
end

cuts = [];
for i = find(p.knapsack.type == 0)
    a = p.knapsack.a{i};
    b = p.knapsack.b{i};
    v = p.knapsack.variables{i};
    x_ = x(v);
    
    cut_ = knapsack_create_cover_cut(a,b,x_,alg,p.gubs);
    if ~isempty(cut_)
        cut = spalloc(1,length(p.c)+1,0);
        cut(1) = cut_(1);
        cut(1 + v) = cut_(2:end);
        cuts = [cuts;cut];
    end
end

