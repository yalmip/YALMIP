function [L,U] = presolve_glover_sherali_cover(row,p)
% Second-order cover inequalities
% Fred Glover · Hanif D. Sherali
% Assumed input format row*[1;x] >=0

% Some fixing of variables has failed miserably
if isinf(row(1))
    L = [];
    U = [];
    return
end

u = p.binarycardinality.up;
d = p.binarycardinality.down;

% Eliminate variables fixed to 0 or 1
k1 = (p.lb(p.binary_variables) == 1);
k2 = (p.ub(p.binary_variables) == 0);
if any(k1) || any(k2)    
    k1 = find(k1);
    row(1) = row(1) + sum(row(1+p.binary_variables(k1)));
    row(1+p.binary_variables(k1)) = 0;
    u = u - length(k1);
    d = d - length(k1);
    k2 = find(k2);    
    row(1+p.binary_variables(k2)) = 0;
    p.binary_variables = setdiff(p.binary_variables,p.binary_variables(union(k1,k2)));    
    if isempty(p.binary_variables)
        L = [];
        U = [];
        return
    end
end

% if length(p.knapsack.a)>0
%     u_ = p.binarycardinality.up;
%     row_hash = sum(p.hash(find(row(2:end))));
%     found = find(row_hash == p.knapsack.hash);
%     for i = 1:length(found)
%         type = p.knapsack.type(found(i));
%         if type >= 2 && type <= 7
%             u_ = min(u_,p.knapsack.b{found(i)});
%         end
%     end
%     [u u_]
%     u = min(u_,u);
% end

% Reduce to binary only. Everything else should have bee fixed
if ~(length(p.c) == length(p.binary_variables))
    row = row([1 1+p.binary_variables]);
end

% This got trivial, all elements can be turned on
% if length(row)-1 <= u
%     if sum(row(2:end)) < -row(1)
%         % Infeasible
%         L = inf(length(p.c),1);
%         U = -inf(length(p.c),1);        
%     else
%         % feasible
%         L = [];
%         U = [];        
%     end
%     return
% end

% Extract knapsack ax>=a0
a = row(2:end);
a0 = -row(1);

% Core method
[L_,U_] = glover_sherali_raw(a,a0,u,d);
    
% Map back to binary positions
L = -inf(length(p.c),1);
U = inf(length(p.c),1);
L(p.binary_variables) = L_;
U(p.binary_variables) = U_;

