function [prob,problem] = yalmip2geometric(options,F_struc,c,Q,K,ub,lb,mt,linear_variables,extended_variables);
   
problem = 0;
prob    = [];
h = [];
A = []; % powers in objective/inequalities
b = []; % coefficients in objective/inequalities
G = []; % powers in equalities
h = []; % coefficients in equalities

% **********************************************************
% Initial sanity check for posynomial problem structure
% **********************************************************
if any(ub<0)
    problem = -8;
    return
end

% **********************************************************
% Wrong quadratic (negative terms)
% Fix: Add code for the case min -x1*x2 ===> min 1/x1x2
% **********************************************************
if any(any(Q<0))
    problem = -8;
    return
end

% **********************************************************
% Setup data related to  objective
% **********************************************************
vars_in_objective = find(c);
A = mt(vars_in_objective,linear_variables);
b = c(vars_in_objective);
map_pos = 0;
map = repmat(map_pos,length(b),1);
 
% **********************************************************
% Invert negative monomial, or this is not a posynomial
% **********************************************************
if any(b<0)
    if nnz(b)==1
        b(find(b)) = -1/b(find(b));
        A = -A;
    else
        problem = -8;
        return
    end
end

% **********************************************************
% Setup data related to inequalities sum(?x^?) > 0
% Loop through all inequalities, find the element with
% positive coefficient, divide by this term. 
% **********************************************************
mte = blkdiag(0,mt); % Extend the monomial table with the monomial x^0
% These are strange, either trivial or not posynomials (multiple terms have
% positive coefficients, such as 1 -x + y > 0 or x + y > 0)
k = find((sum(F_struc(1+K.f:end,:)>0,2))>1);
if ~isempty(k)
    for j = k(:)'
        if all(F_struc(j,:)>=0)
            % Redundant x+y > -something
        elseif all(F_struc(j,:)<=0)
            % Infeasible x+y < -1
            problem = 1;
            return
        else
            % Not posynomial at least
            problem = -8;
            return
        end
    end
end
% remove those redundant constraints
F_struc(K.f + k,:)=[];
k = [];

if any(F_struc(:,1) < 0 & all(F_struc(:,2:end)<=0,2))
    % This is something like x + y < -1
    problem = 1;
    return;
end

% horribly hacked to avoid FIND and sparse subsref inside the loop
% FIXME: Assumes no constraints of the type x + y < 0
ss = full([F_struc > 0]*[(1:size(F_struc,2))]');
[all_non_zeros,temp,vals] = find(F_struc(1+K.f:end,:)');
pos_vals = (vals(vals>0));
nnz_per_row = sum(F_struc~=0,2);
index = 1;
A = A';
map = map';
for j = setdiff(1+K.f:size(F_struc,1),k)   
    k = ss(j);
    vars_in_c = all_non_zeros(index:index + nnz_per_row(j) - 1)';
    Atemp = mte(vars_in_c,linear_variables+1);
    Atemp = Atemp - repmat(mte(k,linear_variables+1),size(Atemp,1),1);
    % f = F_struc(j,:);
    % btemp = reshape(-f(vars_in_c)/f(k),length(vars_in_c),1);
    btemp = reshape(-vals(index:index + nnz_per_row(j) - 1)/pos_vals(1),length(vars_in_c),1);
    % btemp = reshape(-vals(index:index + nnz_per_row(j) - 1)/pos_vals(k),length(vars_in_c),1);
    removed = find(sum(abs(Atemp),2)==0);
    Atemp(removed,:) = [];
    btemp(removed) = [];
    index = index +  nnz_per_row(j);
    pos_vals(1)=[];
    if length(btemp) > 0
        A = [A Atemp'];
        b = [b;full(btemp)];
        map_pos = map_pos + 1;
        map = [map map_pos*ones(1,length(btemp))];
    end
end
map = map';
A = A';

% **********************************************************
% Fix equality constraints coming from fractional powers
% of posynomials.
% An equality constraint a(x) = t can be relaxed to a(x)<t
% if only positive powers of t are used in the program.
% NOTE : YALMIP defines these equalities as t-a(x)==0
% FIX : Is this check enough?
% FIX : Speed things up...
% **********************************************************
for j = 1:1:K.f
    k = find(F_struc(j,:)>0);
    if length(k)>1
        problem = -8;
        return        
    else
        if k==1 | ~ismember(k-1,extended_variables)
            % Monomial equality ok!
        else
            if all(A(:,find(ismember(linear_variables,k-1)))>=0)
            else
                problem = -8;
                return                    
            end
        end
    end
end

% **********************************************************
% Setup data related to inequalities derived from equalities
% i.e. ax^b == 1 replaced with  ax^b<1, (x^-b)/a < 1
% (except for extended variables according to above)
% **********************************************************
for j = 1:1:K.f
    k = find(F_struc(j,:)>0);
    if length(k) == 1
        vars_in_c = find(F_struc(j,:));
        Atemp = mte(vars_in_c,linear_variables+1);
        Atemp = Atemp - repmat(mte(k,linear_variables+1),size(Atemp,1),1);
        btemp = reshape(-F_struc(j,vars_in_c)/F_struc(j,k),length(vars_in_c),1);
        removed = find(sum(abs(Atemp),2)==0);
        Atemp(removed,:) = [];
        btemp(removed) = [];
        
        if ~ismember(k-1,extended_variables)
            if length(btemp)==1
                G = [G;Atemp];
                h = [h;btemp];
            else
                % a(x)+b(x) == c(x) not supported
                problem = -4;
                return
            end
        else
            % Just add upper inequalities
            A = [A;Atemp];
            b = [b;btemp];
            map_pos = map_pos + 1;
            this_map = repmat(map_pos,length(btemp),1);
            map = [map;this_map];
        end
    else
        % a(x) < b(x) + c(x) not supported
        problem = -4;
        return
    end

end

% **********************************************************
% MOSEK does not like upper boud == lower bound
% **********************************************************
if ~(isempty(lb) | isempty(ub))
    fixed_variables = find(lb(linear_variables)==ub(linear_variables));
    if ~isempty(fixed_variables)
        fixed_values = lb(linear_variables(fixed_variables));
        if any(fixed_values==0)
            zeros_vars = find((lb(linear_variables)==ub(linear_variables)) & (lb(linear_variables)==0));
            if any(A(:,zeros_vars)<0)
                problem = 1;
                return
            end
        end
        
        for i = 1:size(A,1)
            this_gain = 1;
            for j = 1:size(fixed_variables)
                this_gain = this_gain*fixed_values(j)^A(i,fixed_variables(j));
                A(i,fixed_variables(j))=0;
            end
            b(i)=b(i)*this_gain;
        end
    end
end
if ~isempty(ub)
    for i = 1:length(ub(linear_variables))
        if ~(isinf(ub(linear_variables(i))) | ub(linear_variables(i))==0)
            A(end+1,i) = 1;
            b(end+1) = 1/ub(linear_variables(i));
            map_pos = map_pos + 1;
            map(end+1)=map_pos;
        end
    end
end
if ~isempty(lb)
    for i = 1:length(lb(linear_variables))
        if ~(isinf(lb(linear_variables(i))) | lb(linear_variables(i))==0)
            A(end+1,i) = -1;
            b(end+1) = lb(linear_variables(i));
            map_pos = map_pos + 1;
            map(end+1)=map_pos;
        end
    end
end

prob.b = b;prob.A = A;prob.map = map;prob.G = G;prob.h = h;
