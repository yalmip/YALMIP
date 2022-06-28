function p = detectAtMost(p)
atmostgroups = {};
atmostbounds = [];
atmostvariables = [];
atleastgroups = {};
atleastbounds = [];
atleastvariables = [];
integers = unique([p.integer_variables p.binary_variables]);
if any(p.K.l) && ~isempty(integers)
    nint = length(integers);
    Aineq = -p.F_struc(p.K.f + (1:p.K.l),2:end);
    bineq = p.F_struc(p.K.f + (1:p.K.l),1);
    notinteger_var_index = setdiff(1:length(p.lb),integers);
    only_integer = ~any(Aineq(:,notinteger_var_index),2);
    Aineq_bin = Aineq(find(only_integer),integers);
    bineq_bin = bineq(find(only_integer),:);
    % Detect groups with constraints #(d_i ~= 0) <= k (for binaries)
    atmostgroups = {};
    for i = 1:size(Aineq_bin,1)
        if bineq_bin(i) >= 1
            [ix,jx,sx] = find(Aineq_bin(i,:));
            % 0/1 or -1/0 variables with with bounded cardinality
            case1 = all(sx < 0) && all(p.ub(integers(jx))<=0);
            case2 = all(sx > 0 )&& all(p.lb(integers(jx))>=0);
            if case1 || case2
                [sorted,loc] = sort(abs(sx),'ascend');
                implied_cardinality = max(find(cumsum(sorted)<=bineq_bin(i)));
                atmostgroups{end+1} = integers(jx);
                atmostbounds = [atmostbounds implied_cardinality];
               % atmostbounds = [atmostbounds floor(bineq_bin(i))];
                atmostvariables = [atmostvariables integers(jx)];
            end        
        elseif bineq_bin(i) <= -1
            [ix,jx,sx] = find(Aineq_bin(i,:));
            % 0/1 or -1/0 variables with with bounded cardinality
            case1 = all(sx == 1) && all(p.lb(integers(jx))==-1) && all(p.ub(integers(jx))==0);
            case2 = all(sx == -1) && all(p.lb(integers(jx))==0) && all(p.ub(integers(jx))==1);               
            if case1 || case2
                atleastgroups{end+1} = integers(jx);
                atleastbounds = [atleastbounds floor(bineq_bin(i))];
                atleastvariables = [atleastvariables integers(jx)];
            end
        end        
    end
end
% Globally detected
p.atmost.groups = atmostgroups;
p.atmost.bounds = atmostbounds;
p.atmost.variables = unique(atmostvariables);
p.atleast.groups = atleastgroups;
p.atleast.bounds = -atleastbounds;
p.atleast.variables = unique(atleastvariables);
% We might use local later
p.localatmost.groups = {};
p.localatmost.bounds = [];
p.localatmost.variables = [];