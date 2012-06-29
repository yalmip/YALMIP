function model = fixbounds(model,these);

if nargin == 1
    polynomials = find((model.variabletype ~= 0));
else
    polynomials = find((model.variabletype ~= 0));
    polynomials = polynomials(find(ismember(polynomials,these)));
end

for i = 1:length(polynomials)
    j = polynomials(i);
    if j<=length(model.lb)
        bound = powerbound(model.lb,model.ub,model.monomtable(j,:));
        model.lb(j) = max(model.lb(j),bound(1));
        model.ub(j) = min(model.ub(j),bound(2));
    end
end
if ~isempty(model.integer_variables)
    model.lb(model.integer_variables) = fix(model.lb(model.integer_variables)-1e-8);
    model.ub(model.integer_variables) = fix(model.ub(model.integer_variables)+1e-8);
end
if ~isempty(model.binary_variables)
    model.lb(model.binary_variables) = fix(model.lb(model.binary_variables)-1e-8);
    model.ub(model.binary_variables) = fix(model.ub(model.binary_variables)+1e-8);
end

