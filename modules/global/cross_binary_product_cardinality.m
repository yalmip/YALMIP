function p = crossBinaryProductandCardinality(p)

if ~isempty(p.binaryProduct)
    for i = find(p.knapsack.type == 2 | p.knapsack.type == 4)
        g1 = ismember(p.binaryProduct(:,2),p.knapsack.variables{i});
        g2 = ismember(p.binaryProduct(:,3),p.knapsack.variables{i});
        s = find(g1 & g2);
        if ~isempty(s)            
            if p.knapsack.b{i}==2
                row = spalloc(1,1+length(p.c),0);
                row(1)=1;
                row(1+p.binaryProduct(s,1))=-1;
                p=addInequality(p,row);
%                p.atmost.groups{end+1}=p.binaryProduct(s,1);
%                p.atmost.bounds(end+1)=1;
%                p.atmost.variables = unique([p.atmost.variables p.binaryProduct(s,1)']);
                
                p.knapsack.a{end+1} = ones(1,length(p.binaryProduct(s,1)));
                p.knapsack.b{end+1} = 1;
                p.knapsack.variables{end+1} = p.binaryProduct(s,1);             
                p.knapsack.type(end+1) = 3;
            end
        end
    end
end