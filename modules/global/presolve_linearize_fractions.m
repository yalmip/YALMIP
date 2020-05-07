function p = presolve_linearize_fractions(p)
% NOT FINISHED YET.
% Handle **1/t** + x/t
% Any elementwise constraints?
if p.K.l >0
    % Any quadratic terms in model?
    if any(p.variabletype==4)
        newA = [];
        newB = [];
        for i = 1:p.K.l
            U = p.F_struc(p.K.f+i,1);
            [aux,col,val] = find(p.F_struc(p.K.f+i,2:end));
            col = col(:);
            % Number of monomials in constraint
            m = length(col);
            % Compute the number of monomials which are divided by a
            % variable xi, and search for the case when all are divided by
            % one variable
            xi = find(sum(full(p.monomtable(col,:)==-1))==m);
            if length(xi)==1 && p.lb(xi) > 0
                % OK, this is something like f1(x)xi + ...fm(x)xi <= q  
                new_var = [];
                for j = 1:m               
                    fj = p.monomtable(col(j),:)
                    fj(xi)=0;
                    new_var = [new_var findrows(p.monomtable,fj)];
                end
                if length(new_var) == m
                end
            end
        end
    end
end