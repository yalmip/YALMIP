function p = detectAtMost(p)
atmostgroups = {};
atmostbounds = [];
atmostvariables = [];
if p.K.l > 0 & ~isempty(p.integer_variables)
    nint = length(p.integer_variables);
    Aineq = -p.F_struc(p.K.f + (1:p.K.l),2:end);
    bineq = p.F_struc(p.K.f + (1:p.K.l),1);
    notinteger_var_index = setdiff(1:length(p.lb),p.integer_variables);
    only_integer = ~any(Aineq(:,notinteger_var_index),2);
    Aineq_bin = Aineq(find(only_integer),p.integer_variables);
    bineq_bin = bineq(find(only_integer),:);
    % Detect groups with constraints #(d_i ~= 0) <= k (for binaries)
    atmostgroups = {};
    for i = 1:size(Aineq_bin,1)
        if bineq_bin(i) >= 1
            [ix,jx,sx] = find(Aineq_bin(i,:));
            % 0/1 or -1/0 variables with with bounded cardinality
            if all(sx == -1) && ((all(p.lb(jx)==-1) && all(p.ub(jx)==0)) || (all(p.ub(jx)==1) && all(p.lb(jx)==0)))                
                atmostgroups{end+1} = p.integer_variables(jx);
                atmostbounds = [atmostbounds floor(bineq_bin(i))];
                atmostvariables = [atmostvariables p.integer_variables(jx)];
            end
        end
    end
end
p.atmost.groups = atmostgroups;
p.atmost.bounds = atmostbounds;
p.atmost.variables = atmostvariables;