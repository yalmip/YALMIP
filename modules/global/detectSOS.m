function p = detectSOS(p)
sosgroups = {};
sosvariables = [];
cardinalitygroups = {};
cardinalityvariables = {};
cardinalitysize = {};
if any(p.K.f) & ~isempty(p.binary_variables)
    nbin = length(p.binary_variables);
    Aeq = -p.F_struc(1:p.K.f,2:end);
    beq = p.F_struc(1:p.K.f,1);
    notbinary_var_index = setdiff(1:length(p.lb),p.binary_variables);
    only_binary = ~any(Aeq(:,notbinary_var_index),2);
    Aeq_bin = Aeq(find(only_binary),p.binary_variables);
    beq_bin = beq(find(only_binary),:);
    % Detect groups with constraints sum(d_i) == 1
    sosgroups = {};
    for i = 1:size(Aeq_bin,1)
        if beq_bin(i) == 1
            [ix,jx,sx] = find(Aeq_bin(i,:));
            if all(sx == 1)
                sosgroups{end+1} = p.binary_variables(jx);
                sosvariables = [sosvariables p.binary_variables(jx)];
            end
        elseif beq_bin(i) > 1
            [ix,jx,sx] = find(Aeq_bin(i,:));
            if all(sx == 1)
                cardinalitygroups{end+1} = p.binary_variables(jx);
                cardinalityvariables = [cardinalityvariables p.binary_variables(jx)];
                cardinalitysize{end+1} = beq_bin(i);
            end
        end
    end
end
p.sosgroups = sosgroups;
p.sosvariables = sosvariables;
p.cardinalitygroups = cardinalitygroups;
p.cardinalityvariables = cardinalityvariables;
p.cardinalitysize = cardinalitysize;

