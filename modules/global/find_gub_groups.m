function [variables, index] = findGUBGroup(p,s)
variables = [];index = [];
for i = find((p.knapsack.type == 3) | (p.knapsack.type == 5))
    if all(ismember(s,p.knapsack.variables{i}))
        variables = p.knapsack.variables{i};
        index = i;
        return
    end
end