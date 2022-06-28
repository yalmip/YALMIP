function covers = knapsack_add_cover_cut(p,x,covers,alg,upper)

if ~isinf(upper) && p.LinearBinaryPositiveCost
    % We can add a global knapsack from c'*x <= upper
    relevant_to_use = find(p.c & p.c <= upper);
    p.knapsack.a{end+1} = p.c(relevant_to_use);
    p.knapsack.b{end+1} = upper;
    p.knapsack.variables{end+1} = relevant_to_use;
    p.knapsack.type(end+1) = 0;    
end

for i = find(p.knapsack.type == 0)
    a = p.knapsack.a{i};
    b = p.knapsack.b{i};
    v = p.knapsack.variables{i};
    x_ = x(v);
    switch alg
        case 'crowder'
            % Heuristics from H. Crowder, E. Johnson, M. Padberg
            % Solving large-scale 0–1 linear programming programs
            [val,loc] = sort((1-x_)./a(:),'ascend');
        case 'gu'
            [val,loc] = sort((x_),'descend');
        otherwise
            error('Unsupport cover cut separation')
    end
    
    % Initial cover
    C = min(find(cumsum(a(loc))>b));
    %  if ~(sum(x_(loc(1:C)))<=(C-1)*1.05)
    
    % Apply Balas lifting
    q = knapsack_cover_lift_balas(a,loc(1:C));
    
    row = spalloc(1,length(p.c)+1,0);
    row(1) = C-1;
    row(1 + v)=-q;
    
    if row*[1;x] < -(C-1)*0.05
        
        covers = [covers;row];
        
        % This is actually a GUB constraint, keep for quick preprocessing
        if row(1)==1
            if all(nonzeros(row(2:end))==-1)
                j = find(row(2:end));
                p.atmost.groups{end+1} = j;
                p.atmost.bounds(end+1) = 1;
                p.atmost.variables = unique([p.atmost.variables find(row(2:end))]);
                p.knapsack.variables{end+1} = j;
                p.knapsack.b{end+1} = 1;
                p.knapsack.a{end+1} = ones(1,length(j));
                p.knapsack.type(end+1) = 5;            
            end
        end
    end
    covers = unique(covers,'rows');
end
