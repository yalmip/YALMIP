function p = propagate_forbidden_region(p)

if ~isempty(p)
    for i = 1:length(p.evalMap)
        if any(p.evalMap{i}.properties.forbidden)
            for spliton = p.evalMap{i}.variableIndex
                if p.ub(spliton) >= p.evalMap{i}.properties.forbidden(1) && p.ub(spliton) <= p.evalMap{i}.properties.forbidden(2)
                    p.ub(spliton) = p.evalMap{i}.properties.forbidden(1);
                    return
                elseif p.lb(spliton) >= p.evalMap{i}.properties.forbidden(1) && p.lb(spliton) <= p.evalMap{i}.properties.forbidden(2)
                    p.lb(spliton) = p.evalMap{i}.properties.forbidden(2);
                    return
                end
            end
        end
    end
end
