function row = fix_continuous_in_row_at_best_possible(row,p,upper)
% Used when dealing with continuous variable in knapsacks
% Just assume it is best possible to get necessary conditions
% Note, p has to be the global model

% A common case is a continuous variable as objective
% If solution available,we can tighten thisbound
if p.UnitContinuousCost && ~isinf(upper)
    k = find(p.c);    
    if p.c(k) == -1
        % Some variable is maximized
        % It will never be lower than current objective
        p.lb(k) = max(p.lb(k),-upper);
    elseif p.c(k) == 1
        % Some variable is minimized
        % It will never be larger than current objective
        p.ub(k) = min(p.ub(k),upper);
    end
end

% We are deriving necessary conditions on binaries
% Fix value of continuous at most favorable and remove
if ~isempty(p.noninteger_variables)
    for i = p.noninteger_variables
        if row(i+1) < 0
            row(1) = row(1) + row(i+1)*p.lb(i);
        elseif row(i+1) > 0
            row(1) = row(1) + row(i+1)*p.ub(i);
        end
    end
    row(1+p.noninteger_variables) = 0;
end
