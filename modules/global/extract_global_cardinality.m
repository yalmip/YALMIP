function p = extract_global_cardinality(p)
% Some code requires a global cardinality constraint
% Look through all knapsacks to find one
% FIXME: Derive an upper bound from groups
p.globalcardinality.up = inf;
p.globalcardinality.down = 0;
for i = 1:length(p.knapsack.variables)
    if p.knapsack.type(i) > 0
        if all(ismember(p.binary_variables,p.knapsack.variables{i}))
            p.globalcardinality.up = min(p.knapsack.b{i},p.globalcardinality.up);
        end
    end
end