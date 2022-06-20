function [p,negated_binary] = detectNegatedBinaries(p)
negated_binary = [];
if ~isempty(p.integer_variables) && nnz(p.Q)==0 && isempty(p.nonlinear)
    negated_binary = find(p.ub(p.integer_variables)==0 & p.lb(p.integer_variables)==-1);
    if ~isempty(negated_binary)
        changed = p.integer_variables(negated_binary);
        p.lb(changed) = 0;
        p.ub(changed) = 1;
        p.c(changed) = -p.c(changed);
        p.F_struc(:,1+changed) = -p.F_struc(:,1+changed);       
        p.binary_variables = union(p.binary_variables,changed);
        p.integer_variables(negated_binary) = [];        
    end
end
