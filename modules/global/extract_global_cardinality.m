function p = extract_global_cardinality(p)
% Derive crude bounds on global cardinality from knapsacks
if length(p.knapsack.b) > 0
    start = [];
    total = 0;
    for i = 1:length(p.knapsack.b)
        rate(i) = p.knapsack.b{i}/sum(p.knapsack.a{i});
    end
    [~,loc] = sort(rate,'ascend');
    A = zeros(0,length(p.c));
    b = [];
    for k = 1:length(p.knapsack.b)
        j = loc(k);       
        if ismember(p.knapsack.type(j),[2 3 4 5])
            A(end+1,p.knapsack.variables{j})=1;
            b(end+1,1)=p.knapsack.b{j};
            s = setdiff(p.knapsack.variables{j},start);
            start = union(start,s);
            total = total + min(p.knapsack.b{j},length(s));
        end
    end
    if length(start)==length(p.binary_variables)
        p.binarycardinality.up = min(p.binarycardinality.up,total);
    else
        p.binarycardinality.up = min(p.binarycardinality.up,total + length(p.binary_variables) - length(start));
    end
    start = [];
    total = 0;
    % Use set covers/lower cardinality
    for k = 1:length(p.knapsack.b)
        if ismember(p.knapsack.type(k),[9 10])
            already_used = intersect(p.knapsack.variables{k},start);
            s = setdiff(p.knapsack.variables{k},start);
            start = union(start,s);
            if length(s)>0
                total = total + min(max(0,p.knapsack.b{k}-length(already_used)));
            end
        end
    end
    %         for j = 1:length(p.atleast.groups)
    %     s = setdiff(p.atleast.groups{j},start);
    %     start = union(start,s);
    %     total = total + min(p.atleast.bounds(j),length(s));
    % end
    p.binarycardinality.down = max(p.binarycardinality.down,total);
    % Silly in case we missed something
    k = nnz(p.lb(p.binary_variables)==1);
    p.binarycardinality.down = max(p.binarycardinality.down,k);
end