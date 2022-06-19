function p = propagate_atmost(p)

for j = 1:length(p.atmost.groups)
    if p.atmost.bounds(j)==1
        xy = p.atmost.groups{j};
        locked = find(p.lb(xy) == 1);
        if length(locked)>0
            if length(locked)>1
                error
            end
            sisters = xy(xy~=xy(locked));
            for k = sisters
                if p.lb(k) > 0
                    error
                    break
                else
                    p.lb(k) = 0;
                    p.ub(k) = 0;
                end
            end
        end
    end
end