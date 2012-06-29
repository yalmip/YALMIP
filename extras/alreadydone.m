function done = alreadydone(variable,method,goal_vexity)
global ALREADY_MODELLED REMOVE_THESE_IN_THE_END
done = 0;
if length(ALREADY_MODELLED) >= variable
    if ~isempty(ALREADY_MODELLED{variable})
        if isequal(ALREADY_MODELLED{variable}.method,method)
            if isequal(ALREADY_MODELLED{variable}.goal_vexity,goal_vexity)
                done = 1;
            end
        elseif isequal(ALREADY_MODELLED{variable}.method,'milp')
            done = 1;
        end
    end
end
