function model = build_recursive_scheme(model);

model.evaluation_scheme = [];
model.deppattern = model.monomtable | model.monomtable;

% Figure out arguments in all polynomials & sigmonials. This info is
% used on several places, so we might just as well save it
model.monomials = find(model.variabletype);
model.monomialMap = cell(length(model.monomials),1);
model.evaluation_scheme = [];
%M = model.monomtable(model.monomials,:);
MM= model.monomtable(model.monomials,:)';
for i = 1:length(model.monomials)
 %   model.monomialMap{i}.variableIndex = find(model.monomtable(model.monomials(i),:));
 %   model.monomialMap{i}.variableIndex = find(M(i,:));
    model.monomialMap{i}.variableIndex = find(MM(:,i));
end

if ~isempty(model.evalMap)

    remainingEvals  = ones(1,length(model.evalVariables));
    remainingMonoms = ones(1,length(model.monomials));
    model = recursive_call(model,remainingEvals,remainingMonoms);

    % Define a dependency structure 
    % (used to speed up Jacobian computations etc)
    model.isevalVariable = zeros(1,size(model.monomtable,1));
    model.isevalVariable(model.evalVariables) = 1;
    compute_depenedency = 1:size(model.monomtable,1);
    % remove purely linear variables, dependency is already computed
    compute_depenedency = setdiff(compute_depenedency,setdiff(find(model.variabletype==0),model.evalVariables));
    if isequal(model.evaluation_scheme{1}.group,'monom')
        % This first level of monomials are simply defined by the monomial
        % table, hence dependency is already computed
        compute_depenedency = setdiff(compute_depenedency,model.monomials(model.evaluation_scheme{1}.variables));        
    end
    for i = compute_depenedency
        k = depends_on(model,i);
        model.deppattern(i,k) = 1;
    end
else
    % Only monomials
    model.evaluation_scheme{1}.group = 'monom';
    model.evaluation_scheme{1}.variables = 1:nnz(model.variabletype);
end

function r = depends_on(model,k)
if model.variabletype(k)
    vars = find(model.monomtable(k,:));
    r = vars;
    for i = 1:length(vars)
        r = [r depends_on(model,vars(i))];
    end
elseif model.isevalVariable(k)%ismember(k,model.evalVariables)
    j = find(k == model.evalVariables);
    r = [];
    for i = 1:length(model.evalMap{j}.variableIndex)
        argument = model.evalMap{j}.variableIndex(i);        
        r = [r argument depends_on(model,argument)];
    end
else
    r = k;
end

function model = recursive_call(model,remainingEvals,remainingMonoms)

if ~any(remainingEvals) & ~any(remainingMonoms)
    return
end

stillE = find(remainingEvals);
stillM = find(remainingMonoms);
    
% Extract arguments in first layer
if any(remainingEvals)
    for i = 1:length(model.evalMap)
        composite_eval_expression(i) = any(ismembcYALMIP(model.evalMap{i}.variableIndex,model.evalVariables(stillE)));
        composite_eval_expression(i) = composite_eval_expression(i) | any(ismembcYALMIP(model.evalMap{i}.variableIndex,model.monomials(stillM)));
    end
end
if any(remainingMonoms)    
    if issorted(model.evalVariables(stillE))
        for i = 1:length(model.monomials)
    %    composite_monom_expression(i) = any(ismember(model.monomialMap{i}.variableIndex,model.monomials(stillM)));
    %    composite_monom_expression(i) = composite_monom_expression(i) | any(ismember(model.monomialMap{i}.variableIndex,model.evalVariables(stillE)));    
        composite_monom_expression(i) = any(ismembcYALMIP(model.monomialMap{i}.variableIndex,model.evalVariables(stillE)));
    end
    else
    for i = 1:length(model.monomials)
    %    composite_monom_expression(i) = any(ismember(model.monomialMap{i}.variableIndex,model.monomials(stillM)));
    %    composite_monom_expression(i) = composite_monom_expression(i) | any(ismember(model.monomialMap{i}.variableIndex,model.evalVariables(stillE)));    
        composite_monom_expression(i) = any(ismember(model.monomialMap{i}.variableIndex,model.evalVariables(stillE)));
    end
    end
end

% Bottom layer
if ~isempty(model.monomials) & any(remainingMonoms)
    if ~isempty(find(~composite_monom_expression & remainingMonoms))
        model.evaluation_scheme{end+1}.group = 'monom';
        model.evaluation_scheme{end}.variables = find(~composite_monom_expression & remainingMonoms);
    end
    remainingMonoms = composite_monom_expression & remainingMonoms;
end

% Bottom layer
if ~isempty(model.evalMap) & any(remainingEvals)
    if ~isempty(find(~composite_eval_expression & remainingEvals));
        model.evaluation_scheme{end+1}.group = 'eval';
        model.evaluation_scheme{end}.variables = find(~composite_eval_expression & remainingEvals);
    end
    remainingEvals = composite_eval_expression & remainingEvals;
end

model = recursive_call(model,remainingEvals,remainingMonoms);