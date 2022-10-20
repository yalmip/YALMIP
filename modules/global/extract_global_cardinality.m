function p = extract_global_cardinality(p)
% Some code requires a global cardinality constraint
% Look through all knapsacks to find one
% FIXME: Derive an upper bound from groups
p.globalcardinality.up = inf;
p.globalcardinality.down = 0;
for i = 1:length(p.knapsack.variables)
    if p.knapsack.type(i) > 0
        if length(p.knapsack.variables{i}) == length(p.binary_variables)
            if all(ismember(p.binary_variables,p.knapsack.variables{i}))
                p.globalcardinality.up = min(p.knapsack.b{i},p.globalcardinality.up);
            end
        end
    end
end
if isinf(p.globalcardinality.up)
    % No simple cardinality knapsack was found
    % Maybe there are sub-groups that we can add up?
    the_sum = 0;
    found = zeros(1,length(p.c));
    for i = 1:length(p.knapsack.variables)
        if p.knapsack.type(i) > 0
            if all(ismember(p.knapsack.variables{i},p.binary_variables))
                if ~any(found(p.knapsack.variables{i}))
                    the_sum = the_sum + p.knapsack.b{i};
                    found(p.knapsack.variables{i}) = 1;
                end
            end
        end
    end   
	p.globalcardinality.up = the_sum + sum(found(p.binary_variables)==0);    
    found = zeros(1,length(p.c));
    for i = length(p.knapsack.variables):-1:1
        if p.knapsack.type(i) > 0
            if all(ismember(p.knapsack.variables{i},p.binary_variables))
                if ~any(found(p.knapsack.variables{i}))
                    the_sum = the_sum + p.knapsack.b{i};
                    found(p.knapsack.variables{i}) = 1;
                end
            end
        end
    end   
	p.globalcardinality.up = min(p.globalcardinality.up,the_sum + sum(found(p.binary_variables)==0));    
end