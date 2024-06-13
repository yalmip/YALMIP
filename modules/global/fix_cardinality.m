function xtemp = fix_cardinality(p,xtemp,x)
would = zeros(1,length(x));
for i = find(ismember(p.knapsack.type,[2 3 4 5]))
    k = p.knapsack.variables{i};
    b = p.knapsack.b{i};
    if nnz(xtemp(k)) > b
        n_should_be_zero = length(k) - b;
        [y,loc] = sort(abs(x(k)),'ascend');
        xtemp(k(loc(1:n_should_be_zero))) = 0;
    end
end
for i = find(ismember(p.knapsack.type,[9 10]))
    k = p.knapsack.variables{i};
    b = p.knapsack.b{i};
    if nnz(xtemp(k)) < b
        n_should_be_one = b;
        [y,loc] = sort(abs(x(k)),'descend');
        xtemp(k(loc(1:n_should_be_one))) = 1;
    end
end