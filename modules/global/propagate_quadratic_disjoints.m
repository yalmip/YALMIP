function p = propagate_quadratic_disjoints(p,spliton)

if any(p.quadraticdisjoints)
    if nargin == 2
        loc = find(spliton == p.quadraticdisjoints(1,:));
        if any(loc)
            p = propagate_one(p,loc);
        end
        return
    end
    for loc = 1:size(p.quadraticdisjoints,2)
        p = propagate_one(p,loc);
    end
end

function p = propagate_one(p,loc)
% x >= u or x <= l
spliton = p.quadraticdisjoints(1,loc);
u = p.quadraticdisjoints(2,loc);
l = p.quadraticdisjoints(3,loc);
if p.lb(spliton) > l
    p.lb(spliton) = max(p.lb(spliton),u);
end
if p.ub(spliton) < u
    p.ub(spliton) = min(p.ub(spliton),l);
end
