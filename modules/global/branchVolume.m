function [V,openVariables] = branchVolume(p,openVariables)

if nargin == 1
    d = p.ub(p.branch_variables)-p.lb(p.branch_variables);
    openVariables = find(d~=0);
    d = d(openVariables);
    if any(d<0)
        % This is infeasible!
        V = 0;
    else
        V = (prod(d))^(1/length(d));
    end
else
    if ~p.feasible
        V = 0;
    else
        d = p.ub(p.branch_variables)-p.lb(p.branch_variables);
        d = d(openVariables); 
        if any(d<0)
            % This has been propagated to infeasibility
            V = 0;
        else
            V = (prod(d))^(1/length(d));
        end
    end        
end