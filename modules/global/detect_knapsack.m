function p = detect_knapsack(p)
binary = p.binary_variables;
p.knapsack.a = {};
p.knapsack.b = {};
p.knapsack.variables = {};
% 0 - general
% 1 - integer
% 2 - Cardinality
% 3 - GUB
% 4 - Derived cardinality from knapsack
% 5 - Derived GUB from knapsack
p.knapsack.type = [];
if any(p.K.l) && ~isempty(binary)
    nint = length(binary);
    Aineq = -p.F_struc(p.K.f + (1:p.K.l),2:end);
    bineq = p.F_struc(p.K.f + (1:p.K.l),1);
    notinteger_var_index = setdiff(1:length(p.lb),binary);
    only_integer = ~any(Aineq(:,notinteger_var_index),2);
    Aineq_bin = Aineq(find(only_integer),binary);
    bineq_bin = bineq(find(only_integer),:);
    % Detect groups with constraints #(d_i ~= 0) <= k (for binaries)    
    for i = 1:size(Aineq_bin,1) 
        if all(Aineq_bin(i,:)>=0)            
            [ix,jx,sx] = find(Aineq_bin(i,:));
            p.knapsack.a{end+1} = sx;
            p.knapsack.b{end+1} = floor(bineq_bin(i));
            p.knapsack.variables{end+1} = binary(jx);
            if all(sx == 1) && bineq_bin(i) == 1
                p.knapsack.type(end+1) = 3;
            elseif all(sx == 1) && bineq_bin(i) > 1
                p.knapsack.type(end+1) = 2;
            else
                p.knapsack.type(end+1) = 0;
            end
            % Add a derived cardinality from general knapsack
            if p.knapsack.type(end) == 0
                sx_ = sort(sx,'ascend');
                C = min(find(cumsum(sx_) > p.knapsack.b{end}));
                p.knapsack.a{end+1} = sx*0+1;
                p.knapsack.b{end+1} = C-1;                
                p.knapsack.variables{end+1} = p.knapsack.variables{end};
                if C-1 == 1
                    p.knapsack.type(end+1) = 5;
                else
                    p.knapsack.type(end+1) = 4;
                end
            end
        end        
    end
end