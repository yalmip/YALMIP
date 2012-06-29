function done = milpalreadydone(i)
global ALREADY_MODELLED REMOVE_THESE_IN_THE_END
done = 0;
if length(ALREADY_MODELLED) >= i
    if ~isempty(ALREADY_MODELLED{i})
        if isequal(ALREADY_MODELLED{i}.method,'milp')
            done = 1;
        end
    end
end
