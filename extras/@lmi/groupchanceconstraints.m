function S = groupchanceconstraints(F)
    S_map = containers.Map('KeyType', 'char', 'ValueType', 'any');
    F = flatten(F);
    
    % Single pass: Identify groups and store indices
    for i = 1:length(F.clauses)
        if ~isempty(F.clauses{i}.confidencelevel)
            g = F.clauses{i}.jointprobabilistic;
            g_key = mat2str(g); % Convert to string for dictionary indexing

            if ~isKey(S_map, g_key)
                S_map(g_key) = []; % Initialize empty array for indices
            end
            S_map(g_key) = [S_map(g_key), i]; % Append index to the group
        end
    end

    % Convert map values to cell array
    keys = S_map.keys;
    S = cell(1, length(keys));
    for i = 1:length(keys)
        s.type = '()';
        s.subs{1} = S_map(keys{i});
        S{i} = subsref(F, s);
    end
end
