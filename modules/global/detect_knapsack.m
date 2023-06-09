function p = detect_knapsack(p)
binary = p.binary_variables;

% Type is saved both in a dedicated knapsack structure with copied data
p.rowtype = -ones(1,p.K.f+p.K.l);
% 0 - general knapsack
% 1 - general integer (not used/detected)
% 2 - Cardinality
% 3 - Set pack/Clique
% 4 - Derived cardinality from knapsack
% 5 - Derived set pack/clique from knapsack
% 6 - Cardinality equality
% 7 - Set partition
% 8 - General equality
% 9 - Set cover
% 10- Cardinality (lower)
if (any(p.K.l) || any(p.K.f)) && ~isempty(binary)
    Aineq = -p.F_struc(1:(p.K.f+p.K.l),2:end);
    bineq = p.F_struc(1:(p.K.f+p.K.l),1);
    notbinary_var_index = setdiff(1:length(p.lb),binary);
    only_binary = ~any(Aineq(:,notbinary_var_index),2);
    binary_rows = find(only_binary);
    Aineq_bin = Aineq(binary_rows,binary);
    bineq_bin = bineq(binary_rows,:);
    % Check every row only containing binary variables
    for i = 1:size(Aineq_bin,1)
        isequality = binary_rows(i)<=p.K.f;
        if isequality && bineq_bin(i)<0
            bineq_bin(i) = -bineq_bin(i);
            Aineq_bin(i,:) = -Aineq_bin(i,:);
            p.F_struc(binary_rows(i),:) = -p.F_struc(binary_rows(i),:);
        end
       
        if all(Aineq_bin(i,:)>=0)
            [ix,jx,sx] = find(Aineq_bin(i,:));
            p.knapsack.a{end+1} = sx;
            p.knapsack.b{end+1} = full(floor(bineq_bin(i)));
            p.knapsack.variables{end+1} = binary(jx);
            p.knapsack.hash(end+1) = sum(p.hash(binary(jx)));
            if all(sx == 1) && bineq_bin(i) == 1
                % sum xi <=(==) 1
                if isequality
                    % Set partition
                    p.knapsack.type(end+1) = 7;
                    p.rowtype(binary_rows(i)) = 7;
                else
                    % Set packing
                    p.knapsack.type(end+1) = 3;
                    p.rowtype(binary_rows(i)) = 3;
                end
                p.cliques.table(binary(jx),end+1)=1;
                p.cliques.hash(end+1) = p.knapsack.hash(end);
                
                % Derive GUB (i.e cardinality on all binaries)
                other = setdiff(p.binary_variables,p.knapsack.variables{end});
                k = 1 + nnz(p.ub(other)==1);
                p.binarycardinality.up = min(k,p.binarycardinality.up);
                
            elseif all(sx == 1) && bineq_bin(i) > 1
                % sum xi <=(==) k                
                % Cardinality constraints
                if isequality
                    p.knapsack.type(end+1) = 6;
                    p.rowtype(binary_rows(i)) = 6;
                    if length(jx) == length(p.binary_variables)
                        % Derive GUB (i.e cardinality on all binaries)
                        p.binarycardinality.up = p.knapsack.b{end};
                        p.binarycardinality.down = p.knapsack.b{end};
                    end
                else
                    p.knapsack.type(end+1) = 2;
                    p.rowtype(binary_rows(i)) = 2;
                    % Derive GUB (i.e cardinality on all binaries)
                    not_used = setdiff(1:length(p.binary_variables), binary(jx));
                    active = nnz(p.ub(not_used)==1);
                    k = p.knapsack.b{end} + active;                 
                    p.binarycardinality.up = min(p.binarycardinality.up,k);
                end
            else
                % sum ai xi <=(==) k
                % General
                if isequality
                    p.knapsack.type(end+1) = 8;
                    p.rowtype(binary_rows(i)) = 8;
                else
                    p.knapsack.type(end+1) = 0;
                    p.rowtype(binary_rows(i)) = 0;
                end
                Cliques = derivecliques(p.knapsack.a{end},p.knapsack.b{end});
                for k = 1:length(Cliques)
                    p.cliques.table(binary(jx(Cliques{k})),end+1) = 1;
                    p.cliques.hash(end+1) = sum(p.hash(binary(jx(Cliques{k}))));
                end
            end
            if p.knapsack.type(end) == 0
                % Add derived knowledge from general knapsack
                % by deriving a cover
                sx_ = sort(sx,'ascend');
                b = p.knapsack.b{end};
                v = p.knapsack.variables{end};
                C = min(find(cumsum(sx_) > b));
                p.knapsack.a{end+1} = sx*0+1;
                p.knapsack.b{end+1} = C-1;
                p.knapsack.variables{end+1} = p.knapsack.variables{end};
                if C-1 == 1
                    % Derived set packing
                    p.knapsack.type(end+1) = 5;
                    p.knapsack.hash(end+1) = sum(p.hash(p.knapsack.variables{end}));
                    p.cliques.table(p.knapsack.variables{end},end+1)=1;
                    p.cliques.hash(end+1) = sum(p.hash(p.knapsack.variables{end}));
                else
                    % Derived cardinality
                    p.knapsack.type(end+1) = 4;
                    p.knapsack.hash(end+1) = sum(p.hash(p.knapsack.variables{end}));
                end
                              
                % Based on this cover, can we update global cardinality?             
                not_used = setdiff(1:length(p.binary_variables), binary(jx));
                active = nnz(p.ub(not_used)==1);
                p.binarycardinality.up = min(p.binarycardinality.up,C-1 + active);                
            end
        elseif all(Aineq_bin(i,:)<=0)
            [ix,jx,sx] = find(Aineq_bin(i,:));
            sx = -sx;
            if all(sx == 1) && bineq_bin(i) == -1
                % Set cover
                p.knapsack.a{end+1} = sx;
                p.knapsack.b{end+1} = full(-bineq_bin(i));
                p.knapsack.variables{end+1} = binary(jx);                
                p.knapsack.type(end+1) = 9;              
                p.knapsack.hash(end+1) = sum(p.hash(binary(jx)));
                p.binarycardinality.down = max(p.binarycardinality.down,1);                                
            elseif all(sx == 1) && bineq_bin(i) < 0
                % Cardinality (lower)
                p.knapsack.a{end+1} = sx;
                p.knapsack.b{end+1} = full(ceil(-bineq_bin(i)));
                p.knapsack.variables{end+1} = binary(jx);
                p.knapsack.type(end+1) = 10;
                p.knapsack.hash(end+1) = sum(p.hash(binary(jx)));
                p.binarycardinality.down = max(p.binarycardinality.down,p.knapsack.b{end});                
            end            
        end
    end
end
p.knapsack.count = p.c'*0;
p.knapsack.contract = ones(1,length(p.c));
for i = 1:length(p.knapsack.type)
    p.knapsack.count(p.knapsack.variables{i})=p.knapsack.count(p.knapsack.variables{i})+1;
    if ismember(p.knapsack.type(i),[2 3 4 5])
        t = p.knapsack.b{i}/length(p.knapsack.variables{i});
        p.knapsack.contract(p.knapsack.variables{i}) = p.knapsack.contract(p.knapsack.variables{i})*t;
    end
end