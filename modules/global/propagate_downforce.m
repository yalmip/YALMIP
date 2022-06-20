function [p,p_feasible] = propagate_downforce(p,p_feasible)
if p_feasible
    for i = 1:length(p.downForce)
        forcing = p.downForce{i}.forcing;
        forced = p.downForce{i}.forced;
        if p.ub(forcing)==0
            if any(p.lb(forced)==1)
                p_feasible = 0;
                return
            else
                p.ub(forced)=0;
            end
        end
    end
end
