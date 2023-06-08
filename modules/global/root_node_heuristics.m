function [x_min,upper] = root_node_heuristics(p,x_min,upper)

% Guess 0, except those fixed to something else
x = zeros(length(p.c),1);
x(p.lb == p.ub) = p.lb(p.lb==p.ub);
x = fix_binary_products(p,x);
upperhere = computecost(p.f,p.c,p.Q,x,p);
if upperhere < upper
    if checkfeasiblefast(p,x,p.options.bnb.feastol)
        x_min = x;
        upper = upperhere;
    else
        [upper,x_min] = upper_from_sdpextension(p,x,upper,x_min);
    end
end

% Guess 1, except those fixed to something else
x = zeros(length(p.c),1);
x(p.binary_variables) = 1;
x(p.lb == p.ub) = p.lb(p.lb==p.ub);
x = fix_binary_products(p,x);
upperhere = computecost(p.f,p.c,p.Q,x,p);
if upperhere < upper
    if checkfeasiblefast(p,x,p.options.bnb.feastol)
        x_min = x;
        upper = upperhere;
    else
        [upper,x_min] = upper_from_sdpextension(p,x,upper,x_min);
    end
end

if ~isempty(p.knapsack.variables)
    % Super greedy, we start by sorting knapsacks
    % So we can pick cheapest from cheapest
    x = zeros(length(p.c),1);
    x(p.lb == p.ub) = p.lb(p.lb==p.ub);
    for i = 1:length(p.knapsack.variables)
        c(i) = min(p.c(p.knapsack.variables{i}));
    end
    [~,cloc] = sort(c,'ascend');
    for i = cloc
        if ismember(p.knapsack.type(i),[2 3 4 5])
            if nnz(x(p.knapsack.variables{i})) == 0
                [~,loc] = sort(p.c(p.knapsack.variables{i}),'ascend');
                x(p.knapsack.variables{i}(loc(1))) = 1;
            end
        end
        x = fix_binary_products(p,x);
        upperhere = computecost(p.f,p.c,p.Q,x,p);
        if upperhere < upper            
            if checkfeasiblefast(p,x,p.options.bnb.feastol)
                x_min = x;
                upper = upperhere;
            else
                [upper,x_min] = upper_from_sdpextension(p,x,upper,x_min);
            end
        end
    end
    % Try expensive first
    x = zeros(length(p.c),1);
    x(p.lb == p.ub) = p.lb(p.lb==p.ub);
    [~,cloc] = sort(c,'descend');
    for i = cloc
        if ismember(p.knapsack.type(i),[2 3 4 5])
            if nnz(x(p.knapsack.variables{i})) == 0
                [~,loc] = sort(p.c(p.knapsack.variables{i}),'descend');
                x(p.knapsack.variables{i}(loc(1))) = 1;
            end
        end
        x = fix_binary_products(p,x);
        upperhere = computecost(p.f,p.c,p.Q,x,p);
        if upperhere < upper            
            if checkfeasiblefast(p,x,p.options.bnb.feastol)
                x_min = x;
                upper = upperhere;
            else
                [upper,x_min] = upper_from_sdpextension(p,x,upper,x_min);
            end            
        end
    end
end