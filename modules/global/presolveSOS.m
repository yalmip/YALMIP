function p = presolveSOS(p)
if ~isempty(p.sosgroups)
    for i = 1:length(p.sosgroups)
        for j = 1:length(p.sosgroups{i})
            % Probe what happens if we activate this binary
            k = p.sosgroups{i}(j);
            other = setdiff(p.sosgroups{i},k);
            if p.ub(k) && ~any(p.lb(p.sosgroups{i}))
                % Still free to be activated
                p_ = p;
                p_.ub(other) = 0;
                p_.lb(k) = 1;
                p_ = smashFixed(p_);
                if ~isempty(fixedInfeasibleEquality(p_))      
                    p.ub(k) = 0;
                end
                % Try de-activate too
                p_ = p;                
                p_.ub(k) = 0;
                p_ = smashFixed(p_);
                if ~isempty(fixedInfeasibleEquality(p_))      
                    p.lb(k) = 1;
                end
            end
        end
    end
    % Prune in case some are presolved completely
    for i = 1:length(p.sosgroups)
        if nnz(p.ub(p.sosgroups{i}))==1
            keep(i) = 0;
        else
            keep(i) = 1;
        end
    end
    if ~all(keep)
        p.sosgroups = {p.sosgroups{find(keep)}};
    end
end