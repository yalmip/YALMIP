function model = update_integer_bounds(model);

if ~isempty(model.integer_variables)
    I = model.integer_variables;
    % Clean up things like 0.00001 <= x <= 0.999999 to 0 <= x <= 1
    % however, don't kill equalities 1 <= x <= 1
    lbfixed = fix(model.lb(I));
    ubfixed = fix(model.ub(I));
    fixlb = find(model.lb(I) == lbfixed);
    fixub = find(model.ub(I) == ubfixed);
    model.lb(I) = ceil(model.lb(I)-1e-4);
    model.ub(I) = floor(model.ub(I)+1e-4);
    if ~isempty(fixlb)
        model.lb(I(fixlb)) = lbfixed(fixlb);
    end
    if ~isempty(fixub)
        model.ub(I(fixub)) = ubfixed(fixub);
    end
    if any(model.lb(I) > model.ub(I))
        model.feasible = 0;
    end
end
if ~isempty(model.binary_variables)
    B = model.binary_variables;
    model.ub(B) = min(model.ub(B),1);
    model.lb(B) = max(model.lb(B),0);
    lbfixed = fix(model.lb(B));
    ubfixed = fix(model.ub(B));
    fixlb = find(model.lb(B) == lbfixed);
    fixub = find(model.ub(B) == ubfixed);
    model.lb(B) = ceil(model.lb(B)-1e-4);
    model.ub(B) = floor(model.ub(B)+1e-4);
    if ~isempty(fixlb)
        model.lb(B(fixlb)) = lbfixed(fixlb);
    end
    if ~isempty(fixub)
        model.ub(B(fixub)) = ubfixed(fixub);
    end
    if any(model.lb(B) > model.ub(B))
        model.feasible = 0;
    end
end