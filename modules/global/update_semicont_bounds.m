function p = update_semicont_bounds(p)
if ~isempty(p.semicont_variables)
    redundant = find(p.lb<=0 & p.ub>=0);
    p.semicont_variables = setdiff(p.semicont_variables,redundant);
    % Now relax the model and generate hull including 0
    p.semibounds.lb = p.lb(p.semicont_variables);
    p.semibounds.ub = p.ub(p.semicont_variables);
    p.lb(p.semicont_variables) = min(p.lb(p.semicont_variables),0);
    p.ub(p.semicont_variables) = max(p.ub(p.semicont_variables),0);
end
