function p = detectAtMost(p)
atmostgroups = {};
atmostbounds = [];
atmostvariables = [];
integers = unique([p.integer_variables p.binary_variables]);
if any(p.K.l) & ~isempty(integers)
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
            case1 = all(sx == -1) && all(p.lb(integers(jx))==-1) && all(p.ub(integers(jx))==0);
            case2 = all(sx ==  1) && all(p.lb(integers(jx))==0) && all(p.ub(integers(jx))==1);               
            if case1 || case2
                atmostgroups{end+1} = integers(jx);
                atmostbounds = [atmostbounds floor(bineq_bin(i))];
                atmostvariables = [atmostvariables integers(jx)];
            end
        end
    end
end
p.atmost.groups = atmostgroups;
p.atmost.bounds = atmostbounds;
p.atmost.variables = atmostvariables;
%p = presolve_AtmostFromBinaryproduct(p)

function p = presolve_AtmostFromBinaryproduct(p)
% presolve y = x1*x2 where x1 x2 binary
% thus y is binary
ff = p.F_struc(1+p.K.f:p.K.f+p.K.l,:);
% Search for y >= x1 + x2 - 1
candidates = find(sum(abs(ff),2) == 4 & sum(ff,2) == 0 & ff(:,1)==1);
lower_and = {};
for i = 1:length(candidates)
    row = p.F_struc(candidates(i),2:end);
    xx = find(row==-1);
    yy = find(row==1);
    if length(xx)==2 && length(yy)==1
        if all(ismember(xx,p.binary_variables))
            lower_and{end+1}.x = xx;
            lower_and{end}.y = yy;
        end
    end
end
for i = 1:length(lower_and)
    x = lower_and{i}.x;
    y = lower_and{i}.y;
    y_bounded_by_x1 = 0;
    q = find(ff(:,1) == 0 & ff(:,x(1)+1)==1);
    for j = 1:length(q)
        s = ff(q(j),2:end);
        s(x(1))=0;
        % y + sum z <= x1   (z non-negative)
        if all(s<=0) & s(y)==-1 & all(p.lb(find(s))>=0)
          p.atmost.groups{end+1} = find(s);
          p.atmost.bounds(end+1) = 1;          
        end
    end
    q = find(ff(:,1) == 0 & ff(:,x(2)+1)==1);
    for j = 1:length(q)
        s = ff(q(j),2:end);
        s(x(2))=0;
        if all(s<=0) && s(y) && all(p.lb(find(s)) >=0)
          p.atmost.groups{end+1} = find(s);
          p.atmost.bounds(end+1) = 1;   
        end
    end   
end
