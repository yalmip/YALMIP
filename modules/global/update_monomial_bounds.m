function model = update_monomial_bounds(model,these);

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


