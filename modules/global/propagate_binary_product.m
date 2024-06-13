function p = propagate_binary_product(p)

if isempty(p.binaryProduct)
    return
end
if p.feasible
    % y = x1 & x2 detected in presolve so tighten
    % Any obvious inconsistency?
    x1_true = p.lb(p.binaryProduct(:,2))==1;
    x2_true = p.lb(p.binaryProduct(:,3))==1;
    x1_false = p.ub(p.binaryProduct(:,2))==0;
    x2_false = p.ub(p.binaryProduct(:,3))==0;
    y_false = p.ub(p.binaryProduct(:,1))==0;
    y_true = p.lb(p.binaryProduct(:,1))==1;
    
    if any(x1_true & x2_true & y_false)
        p.feasible = 0;
        return
    end
    if any((x1_false | x2_false) & y_true)
        p.feasible = 0;
        return
    end
    
    if any(x1_true & x2_true & ~y_true) || any(y_true & ~(x1_true & x2_true))
        for i = 1:length(p.binaryProduct)
            % Convoluted code to avoid rare reference
            if x1_true(i)
                if x2_true(i)
                    y = p.binaryProduct(i,1);
                    p.lb(y) = 1;
                end
            elseif (x1_false(i) || x2_false(i)) & ~y_false(i)
                y = p.binaryProduct(i,1);
                p.ub(y) = 0;
            end
            if y_true(i)
                y = p.binaryProduct(i,1);
                x1 = p.binaryProduct(i,2);
                x2 = p.binaryProduct(i,3);
                p.lb(x1) = 1;
                p.lb(x2) = 1;
            end
        end
    end
end