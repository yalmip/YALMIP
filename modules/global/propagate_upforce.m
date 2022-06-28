function p = propagate_upforce(p)

% Find inconsistency of y <= x1 + x2 + ...
if p.feasible
    for i = 1:length(p.upForce)
        forcing = p.upForce{i}.forcing;       
        if p.lb(forcing)==1
            forced = p.upForce{i}.forced;
            if all(p.ub(forced)==0)
                p.feasible = 0;
                return
            end
        end
    end
end
