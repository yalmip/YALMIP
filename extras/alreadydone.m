function done = alreadydone(variable,method,goal_vexity)
global ALREADY_MODELLED ALREADY_MODELLED_INDEX REMOVE_THESE_IN_THE_END
done = 0;
index = find(ALREADY_MODELLED_INDEX == variable);
if ~isempty(index)
    if ~isempty(ALREADY_MODELLED{index})
        if isequal(ALREADY_MODELLED{index}.method,method)
            if isequal(ALREADY_MODELLED{index}.goal_vexity,goal_vexity)
                done = 1;
            end
        elseif isequal(ALREADY_MODELLED{index}.method,'milp') | isequal(ALREADY_MODELLED{index}.method,'integer') | isequal(ALREADY_MODELLED{index}.method,'exact')
            done = 1;
        end
    end
end
% if length(ALREADY_MODELLED) >= variable
%     if ~isempty(ALREADY_MODELLED{variable})
%         if isequal(ALREADY_MODELLED{variable}.method,method)
%             if isequal(ALREADY_MODELLED{variable}.goal_vexity,goal_vexity)
%                 done = 1;
%             end
%         elseif isequal(ALREADY_MODELLED{variable}.method,'milp')
%             done = 1;
%         end
%     end
% end
