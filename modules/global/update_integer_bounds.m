function model = update_integer_bounds(model);

if ~isempty(model.integer_variables)
    model.lb(model.integer_variables) = fix(model.lb(model.integer_variables)-1e-4);
    model.ub(model.integer_variables) = fix(model.ub(model.integer_variables)+1e-4);
end
if ~isempty(model.binary_variables)
    model.lb(model.binary_variables) = fix(model.lb(model.binary_variables)-1e-4);
    model.ub(model.binary_variables) = fix(model.ub(model.binary_variables)+1e-4);
end
if any(model.lb(model.binary_variables) > model.ub(model.binary_variables))
    model.feasible = 0;
end
if any(model.lb(model.integer_variables) > model.ub(model.integer_variables))
    model.feasible = 0;
end