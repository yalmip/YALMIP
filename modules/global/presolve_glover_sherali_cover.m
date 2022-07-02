function [L,U] = presolve_glover_sherali_cover(row,p)
% Second-order cover inequalities
% Fred Glover · Hanif D. Sherali

% Some fixing of variables has failed miserably
if isinf(row(1))
    L = [];
    U = [];
    return
end

u = p.globalcardinality.up;

% Eliminate variables fixed to 1
k1 = (p.lb(p.binary_variables) == 1);
k2 = (p.ub(p.binary_variables) == 0);
if any(k1) || any(k2)    
    k1 = find(k1);
    row(1) = row(1) + sum(row(1+p.binary_variables(k1)));
    row(1+p.binary_variables(k1)) = 0;
    u = u - length(k1);
    k2 = find(k2);
    row(1+p.binary_variables(k2)) = 0;
    p.binary_variables = setdiff(p.binary_variables,p.binary_variables(union(k1,k2)));    
    if isempty(p.binary_variables)
        L = [];
        U = [];
        return
    end
end

% Reduce to binary only. Everything else should have bee fixed
if ~(length(p.c) == length(p.binary_variables))
    row = row([1 1+p.binary_variables]);
end

% This got trivial, all elements can be turned on
if length(row)-1 == u
    if sum(row(2:end)) < -row(1)
        % Infeasible
        L = inf(length(p.c),1);
        U = -inf(length(p.c),1);        
    else
        % feasible
        L = [];
        U = [];        
    end
    return
end

% Extract knapsack ax>=a0
a = row(2:end);
a0 = -row(1);
    
L_ = -inf(length(a),1);
U_ = inf(length(a),1);
% Map to sorted
[a_, loc] = sort(a,'descend');
n = length(a_);

%From here, we only work with binary and in sorted
for j = 1:n
    s = a_(1:j);
    SN(j) = sum(s(s>0));
end
%SN = cumsum(max(a_,0));
try
    SN_u_plus_1 = sum(a_(1:u+1));
catch
    1
end
SN_u_minus_1 = sum(a_(1:u-1));

fixed_at_one = find(a_ > SN_u_plus_1-a0);
if ~isempty(fixed_at_one)
    L_(fixed_at_one) = 1;
    U_(fixed_at_one) = 1;
end

jhat = u+1:n;
fixed_at_zero = jhat(find(a_(jhat) < a0-SN_u_minus_1));
if ~isempty(fixed_at_zero)
    U_(fixed_at_zero) = 0;
    L_(fixed_at_zero) = 0;
end

% Map back to unsorted, and all variables
L = -inf(length(p.c),1);
U = inf(length(p.c),1);
try
    L(p.binary_variables(loc)) = L_;
    U(p.binary_variables(loc)) = U_;
catch
    1
end
