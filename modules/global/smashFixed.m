function p = smashFixed(p,action)
if ~isempty(p.F_struc)
    r = find(p.lb == p.ub);
    if ~isempty(r)
        p.F_struc(:,1) = p.F_struc(:,1) + p.F_struc(:,1+r)*p.lb(r);
        if nargin == 1 || strcmp(action,'keep')
            p.F_struc(:,1+r) = 0;
        else
            p.F_struc(:,1+r) = [];
        end
    end
end