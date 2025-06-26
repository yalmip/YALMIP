function p = presolve_cliques(p)

if size(p.cliques,2)>0 && p.K.l>0 && p.K.f == 0 && nnz(p.K.s)==0 && nnz(p.Q) == 0 && p.feasible
    
    for clique = 1:size(p.cliques.table,2)
        var = find(p.cliques.table(:,clique));
        A = -p.F_struc(:,2:end);
        for i = 1:length(var)-1
            var1 = var(i);
            for j = i+1:length(var)
                var2 = var(j);
                % Should/can var be fixed to 0
                % It's more expensive, and restricting, so since
                % only one can be 1, don't use this one
                var1_expensive = p.c(var1) > p.c(var2);
                var2_expensive = p.c(var2) > p.c(var1);
                equal_expensive = p.c(var2) == p.c(var1);
                var1_restricting = all(A(:,var1) >= A(:,var2));
                var2_restricting = all(A(:,var2) >= A(:,var1));
                if (var1_expensive || equal_expensive) && var1_restricting && p.lb(var1) == 0
                    p.ub(var1) = 0;
                elseif (var2_expensive || equal_expensive) && var2_restricting && p.lb(var2) == 0
                    p.ub(var2) = 0;
                end
            end
        end
    end
end