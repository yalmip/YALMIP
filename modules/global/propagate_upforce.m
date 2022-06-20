function [p,p_feasible] = propagate_upforce(p,p_feasible)

for i = 1:length(p.upForce)
    forcing = p.upForce{i}.forcing;
    forced = p.upForce{i}.forced;
    if p.lb(forcing)==1
        if all(p.ub(forced)==0)
            p_feasible = 0;
            return
        end
    end
end
    