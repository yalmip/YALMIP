function model = presolve_bounds_from_domains(model);

% Sigmonial with non-integer powers must be positive
sigmonials = find((model.variabletype == 4));
for i = 1:length(sigmonials)
    j = sigmonials(i);
    involved = find(model.monomtable(j,:));
    fractional = involved(model.monomtable(j,involved)~=fix(model.monomtable(j,involved)));
    for k = 1:length(fractional)
        model.lb(fractional(k)) = max([1e-9  model.lb(fractional(k))]);
    end
end

% The evaluation based operators can communicate the domain (although
% this really should be available via domain constraints anyway)
% The operator model can however also include range constraints, which we
% just as well might extract and add to the model
for i = 1:length(model.evalVariables)
    j = model.evalVariables(i);
    model.lb(j) = max([model.lb(j) model.evalMap{i}.properties.range(1)]);
    model.ub(j) = min([model.ub(j) model.evalMap{i}.properties.range(2)]);
    j = model.evalMap{i}.variableIndex;
    model.lb(j) = max([model.lb(j) repmat(model.evalMap{i}.properties.domain(1),length(j),1)],[],2);
    model.ub(j) = min([model.ub(j) repmat(model.evalMap{i}.properties.domain(2),length(j),1)],[],2);
end


