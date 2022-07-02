function [p,negated_binary] = detectNegatedBinaries(p)
negated_binary = [];
if ~isempty(p.integer_variables) && nnz(p.Q)==0 && isempty(p.nonlinear)
    negated_binary = find(p.ub(p.integer_variables)==0 & p.lb(p.integer_variables)==-1);
    if ~isempty(negated_binary)
        negated_binary = p.integer_variables(negated_binary);
        p.ub(negated_binary) = min(1,-p.lb(negated_binary));
        p.lb(negated_binary) = 0;        
        p.c(negated_binary) = -p.c(negated_binary);
        p.F_struc(:,1+negated_binary) = -p.F_struc(:,1+negated_binary);       
        p.binary_variables = union(p.binary_variables,negated_binary);
        p.integer_variables = setdiff(p.integer_variables,p.binary_variables);
        if ~isempty(p.x0)
            p.x0(negated_binary) = -p.x0(negated_binary);
        end
    end
end
