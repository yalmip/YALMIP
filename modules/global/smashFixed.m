function p = smashFixed(p)
if ~isempty(p.F_struc) 
    r = find(p.lb == p.ub);
    if ~isempty(r)
        p.F_struc(:,1) = p.F_struc(:,1) + p.F_struc(:,1+r)*p.lb(r);
        p.F_struc(:,1+r) = 0;   
    end
end