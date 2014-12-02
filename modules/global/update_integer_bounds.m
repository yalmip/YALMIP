function model = update_integer_bounds(model);

if ~isempty(model.integer_variables)
    % Clean up things like 0.00001 <= x <= 0.999999 to 0 <= x <= 1
    % however, don't kill equalities 1 <= x <= 1
    lbfixed = fix(model.lb(model.integer_variables));
    ubfixed = fix(model.ub(model.integer_variables));
    fixlb = find(model.lb(model.integer_variables) == lbfixed);
    fixub = find(model.ub(model.integer_variables) == ubfixed);
    model.lb(model.integer_variables) = ceil(model.lb(model.integer_variables)-1e-4);
    model.ub(model.integer_variables) = floor(model.ub(model.integer_variables)+1e-4);
    if ~isempty(fixlb)
        model.lb(model.integer_variables(fixlb)) = lbfixed(fixlb);
    end
    if ~isempty(fixub)
        model.ub(model.integer_variables(fixub)) = ubfixed(fixub);
    end
end
if ~isempty(model.binary_variables)
    lbfixed = fix(model.lb(model.binary_variables));
    ubfixed = fix(model.ub(model.binary_variables));
    fixlb = find(model.lb(model.binary_variables) == lbfixed);
    fixub = find(model.ub(model.binary_variables) == ubfixed);
    model.lb(model.binary_variables) = ceil(model.lb(model.binary_variables)-1e-4);
    model.ub(model.binary_variables) = floor(model.ub(model.binary_variables)+1e-4);
    if ~isempty(fixlb)
        model.lb(model.binary_variables(fixlb)) = lbfixed(fixlb);
    end
    if ~isempty(fixub)
        model.ub(model.binary_variables(fixub)) = ubfixed(fixub);
    end
end
if any(model.lb(model.binary_variables) > model.ub(model.binary_variables))
    model.feasible = 0;
end
if any(model.lb(model.integer_variables) > model.ub(model.integer_variables))
    model.feasible = 0;
end