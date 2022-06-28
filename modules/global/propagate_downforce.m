function p = propagate_downforce(p)

% Find inconsistency of x1 + x2 + ... <= y
if p.feasible
    for i = 1:length(p.downForce)
        forcing = p.downForce{i}.forcing;        
        if p.ub(forcing)==0
            forced = p.downForce{i}.forced;
            if any(p.lb(forced)==1)
                p.feasible = 0;
                return
            else
                p.ub(forced)=0;
            end
        end
    end
end
