% *************************************************************************
% Code for setting the numerical values of nonlinear terms
% *************************************************************************
function p = updateonenonlinearbound(p,changed_var)
if ~isempty(p.bilinears)
    impactedVariables = find((p.bilinears(:,2) == changed_var) | (p.bilinears(:,3) == changed_var));
    x = p.bilinears(impactedVariables,2);
    y = p.bilinears(impactedVariables,3);
    z = p.bilinears(impactedVariables,1);
    x_lb = p.lb(x);
    x_ub = p.ub(x);
    y_lb = p.lb(y);
    y_ub = p.ub(y);
    bounds = [x_lb.*y_lb x_lb.*y_ub x_ub.*y_lb x_ub.*y_ub];
    p.lb(z) = max([p.lb(z) min(bounds,[],2)],[],2);
    p.ub(z) = min([p.ub(z) max(bounds,[],2)],[],2)';
    p.lb(impactedVariables(x==y)<0) = 0;
end
