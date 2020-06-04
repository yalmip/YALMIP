function p = updatemonomialbounds(p);
LU = [p.lb p.ub];
if ~isempty(p.bilinears)
    x = p.bilinears(:,2);
    y = p.bilinears(:,3);
    z = p.bilinears(:,1);
    x_lb = p.lb(x);
    x_ub = p.ub(x);
    y_lb = p.lb(y);
    y_ub = p.ub(y);
    bounds = [x_lb.*y_lb x_lb.*y_ub x_ub.*y_lb x_ub.*y_ub];
    new_lb = max([p.lb(z) min(bounds,[],2)],[],2);
    new_ub = min([p.ub(z) max(bounds,[],2)],[],2);
    % Avoid updating small bounds (numerical reasons)
    update = find(p.lb(z) < p.ub(z)-1e-4);
    p.lb(z(update)) = new_lb(update);
    p.ub(z(update)) = new_ub(update);
    
    implied_pos = find(p.lb(z)>0 & p.lb(x)>=0 & p.lb(y)<0);
    if ~isempty(implied_pos)
        p.lb(y(implied_pos))=0;
    end
    implied_pos = find(p.lb(z)>0 & p.lb(y)>=0 & p.lb(x)<0);
    if ~isempty(implied_pos)
        p.lb(x(implied_pos))=0;
    end
    
    from_square = find(x == y & ~isinf(p.ub(z)));
    if ~isempty(from_square)
        p.ub(x(from_square)) = min(p.ub(x(from_square)),sqrt(p.ub(z(from_square))));        
        p.lb(x(from_square)) = max(p.lb(x(from_square)),-sqrt(p.ub(z(from_square))));        
    end
    
    p = update_integer_bounds(p);
    
    quadratic_variables = p.bilinears(x==y,1);
    p.lb(quadratic_variables(p.lb(quadratic_variables)<0)) = 0;
end

if ~isempty(p.high_monom_model)
    for i = 1:size(p.high_monom_model,1)
        j = p.high_monom_model(i,1);
        % x(j) = x(b)^c
        [a,b,c] = find(p.high_monom_model(i,2:end));
        if length(a) == 1
            if  even(c)
                p.ub(b) = min([p.ub(b) p.ub(j)^(1/c)]);
                p.lb(b) = max([p.lb(b) -p.ub(j)^(1/c)]);
                if p.lb(b) > 0
                    p.lb(j) = max([p.lb(j) p.lb(b)^c]);
                    p.lb(b) = max([p.lb(b) p.lb(j)^(1/c)]);
                elseif p.ub(b) < 0
                    p.lb(b) = max([p.lb(b) -p.ub(j)^(1/c)]);
                    p.ub(b) = min([p.ub(b) -p.lb(j)^(1/c)]);
                end
            else
                if c>0
                    p.ub(b) = min([p.ub(b) sign(p.ub(j))*(abs(p.ub(j))^(1/c))]);
                    p.lb(b) = max([p.lb(b) sign(p.lb(j))*(abs(p.lb(j))^(1/c))]);
                else
                end
            end
        end
    end
end

if ~isequal(LU,[p.lb p.ub])
    p.changedbounds = 1;
end
