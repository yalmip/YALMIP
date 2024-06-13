function p = cross_binary_product_cardinality(p)
% If we have a bunch of yk = xi*xj with x1+x2+x3+x4 + ... <= 2
% it follows that at most 1 yk can be non-zero
if ~isempty(p.binaryProduct)
    for i = find(p.knapsack.type == 2 | p.knapsack.type == 3 | p.knapsack.type == 4)
        g1 = ismember(p.binaryProduct(:,2),p.knapsack.variables{i});
        g2 = ismember(p.binaryProduct(:,3),p.knapsack.variables{i});
        s = find(g1 & g2);
        if ~isempty(s)
            if p.knapsack.b{i}==2
                row = spalloc(1,1+length(p.c),0);
                row(1) = 1;
                row(1+p.binaryProduct(s,1))=-1;
                p=addInequality(p,row);
                p.knapsack.a{end+1} = ones(1,length(p.binaryProduct(s,1)));
                p.knapsack.b{end+1} = 1;
                p.knapsack.variables{end+1} = p.binaryProduct(s,1);
                p.knapsack.type(end+1) = 3;
                p.cliques.table(p.knapsack.variables{end},end+1) = 1;
                p.cliques.hash(end+1) = sum(p.hash(p.knapsack.variables{end}));
            elseif p.knapsack.b{i} == 1
                p.ub(p.binaryProduct(s,1)) = 0;
            end
        end
    end
end