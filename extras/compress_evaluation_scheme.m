function model = compress_evaluation_scheme(model);
scalars = {'exp','log','sin','cos','log2','log10','slog','inverse_internal2'};
for i = 1:length(model.evaluation_scheme)
    if strcmp(model.evaluation_scheme{i}.group,'eval')
        clear fun
        for k = 1:length(scalars)
            for j = 1:length(model.evaluation_scheme{i}.variables)
                fun{k}(j) = strcmp(model.evalMap{model.evaluation_scheme{i}.variables(j)}.fcn,scalars{k});
            end
        end
        for k = 1:length(scalars)
            fun_i = find(fun{k});
            if length(fun_i) > 1
                all_outputs = [];
                all_inputs  = [];
                for j = fun_i
                    all_outputs = [all_outputs model.evalMap{model.evaluation_scheme{i}.variables(j)}.computes];
                    all_inputs  = [all_inputs model.evalMap{model.evaluation_scheme{i}.variables(j)}.variableIndex];
                end
                model.evalMap{model.evaluation_scheme{i}.variables(fun_i(1))}.computes = all_outputs;
                model.evalMap{model.evaluation_scheme{i}.variables(fun_i(1))}.variableIndex = all_inputs;
                model.evaluation_scheme{i}.variables(fun_i(2:end)) = nan;
            end
        end
        model.evaluation_scheme{i}.variables(isnan(model.evaluation_scheme{i}.variables)) = [];    
                
        [model,vars] = removeDuplicates(model,model.evaluation_scheme{i}.variables);
        model.evaluation_scheme{i}.variables = vars;
    end
end

function [model,variables] = removeDuplicates(model,variables)
% optimizer_ evaluation objects leads to multiple copies of similiar
% evaluations which can be performed in one shot
remove = zeros(1,length(model.evalMap));
for i = 1:length(variables)
    for j = i+1:length(variables)
        if ~remove(j)
            if isequal(model.evalMap{variables(i)}.computes,model.evalMap{variables(j)}.computes)
                if  isequal(model.evalMap{variables(i)}.variableIndex,model.evalMap{variables(j)}.variableIndex)
                    remove(j) = 1;
                end
            end
        end
    end
end
if any(remove)
    model.evalMap = {model.evalMap{find(~remove)}};
    variables = variables(find(~remove));
end
