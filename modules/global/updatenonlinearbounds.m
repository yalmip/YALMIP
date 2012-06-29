function p = updatenonlinearbounds(p,changed_var,keepbest);

dbstack
error
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

    p.lb(p.integer_variables) = fix(p.lb(p.integer_variables));
    p.ub(p.integer_variables) = fix(p.ub(p.integer_variables));
    p.lb(p.binary_variables) = fix(p.lb(p.binary_variables));
    p.ub(p.binary_variables) = fix(p.ub(p.binary_variables));

    quadratic_variables = p.bilinears(x==y,1);
    p.lb(quadratic_variables(p.lb(quadratic_variables)<0)) = 0;
end

if ~isempty(p.high_monom_model)
    for i = 1:size(p.high_monom_model,1)
        j = p.high_monom_model(i,1);
        [a,b,c] = find(p.high_monom_model(i,2:end));
        if length(a) == 1
            if even(c)
                % fix this case...
                %   p.ub(b) = min([p.ub(b) sign(p.ub(j))*(abs(p.ub(j))^(1/c))]);
                %   p.lb(b) = max([p.lb(b) sign(p.lb(j))*(abs(p.lb(j))^(1/c))]);
            else
                if c>0
                    p.ub(b) = min([p.ub(b) sign(p.ub(j))*(abs(p.ub(j))^(1/c))]);
                    p.lb(b) = max([p.lb(b) sign(p.lb(j))*(abs(p.lb(j))^(1/c))]);
                else
                    %                     if p.lb(b)>0
                    %                         p.ub(b) = min([p.ub(b) sign(p.lb(j))*(abs(p.lb(j))^(1/c))]);
                    %                         p.lb(b) = max([p.lb(b) sign(p.ub(j))*(abs(p.ub(j))^(1/c))]);
                    %                     end
                end
            end
        end
    end
end