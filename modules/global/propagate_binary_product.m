function [p,p_feasible] = propagate_binary_product(p)
p_feasible = 1;
if isempty(p.binaryProduct)
    return
end
% y = x1 & x2 detected in presolve so tighten
for i = 1:length(p.binaryProduct)
    y = p.binaryProduct(i,1);
    x1 = p.binaryProduct(i,2);
    x2 = p.binaryProduct(i,3);
    % x1 and x2 are true, hence y is true
    if p.lb(x1)==1 && p.lb(x2)==1
        if p.ub(y)==0
            p_feasible = 0;
            return
        else
            p.lb(y) = 1;
        end
    end
    % x1 or x2 are false, hence y is false
    if p.ub(x1)==0 || p.ub(x2)==0
        if p.lb(y)==1
            p_feasible = 0;
            return
        else
            p.ub(y) = 0;
        end
    end
    % y is true, hence x1 and x2 are tru
    if p.lb(y) == 1
        if p.ub(x1)==0 || p.ub(x2)==0
            p_feasible = 0;
            return
        else
            p.lb(x1) = 1;
            p.lb(x2) = 1;
        end
    end
end