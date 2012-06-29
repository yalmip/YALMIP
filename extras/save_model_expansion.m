function done = save_model_expansion(method,goal_vexity,F_graph,properties)
global ALREADY_MODELLED REMOVE_THESE_IN_THE_END ALREADY_MODELLED_INDEX

done = 0;

for variables = properties.models
    save_place = find(variables == ALREADY_MODELLED_INDEX);
    if ~isempty(save_place)%length(ALREADY_MODELLED)>=variables & ~isempty(ALREADY_MODELLED{variables})
        % Ouch, we have already modelled this-one
        if strcmpi(ALREADY_MODELLED{save_place}.method,method)
            % Ok, already modelled using same approach
            done = 1;
            return
        elseif strcmpi(ALREADY_MODELLED{save_place}.method,'graph') & (strcmpi(method,'callback') | strcmpi(method,'integer') | strcmpi(method,'exact'))
            % Replace old graph model with exact model
            REMOVE_THESE_IN_THE_END = [REMOVE_THESE_IN_THE_END  ALREADY_MODELLED{save_place}.index];
            ALREADY_MODELLED{save_place}.goal_vexity = goal_vexity;
            ALREADY_MODELLED{save_place}.method = method;
            ALREADY_MODELLED{save_place}.index  = getlmiid(F_graph);
            ALREADY_MODELLED{save_place}.properties  = properties;
        elseif ((strcmpi(ALREADY_MODELLED{save_place}.method,'integer') | strcmpi(ALREADY_MODELLED{save_place}.method,'exact') | strcmpi(ALREADY_MODELLED{save_place}.method,'callback')) & strcmpi(method,'graph'))
            % Keep old stuff, we are done
            done = 1;
            return
        end
    else
        ALREADY_MODELLED{end+1}.goal_vexity = goal_vexity;
        ALREADY_MODELLED{end}.method = method;
        ALREADY_MODELLED{end}.index  = getlmiid(F_graph);
        ALREADY_MODELLED{end}.properties = properties;
        ALREADY_MODELLED_INDEX = [ALREADY_MODELLED_INDEX;variables];
    end
end
