function p = propagate_bounds_from_evaluations(p);

if ~isempty(p.evalVariables)
    for i = 1:length(p.evalMap)
        p = update_one_eval_bound(p,i);
        p = update_one_inverseeval_bound(p,i);
    end
end