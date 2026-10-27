function upper_bound = upper_from_sdpcuts(p,sdpCuts,upper)
% In models here we have a simple unit continuous objective
% we can use sdp cuts collected on the way to quickly derive
% a crude lower bound
% If this bound is larger than the incumbent upper bound, there is no
% way the optimal solution is here, so we can give up

% Fix given variables
S = sdpCuts.F_struc;
k = find(p.ub == p.lb);
s = S(:,1+k)*p.ub(k);
S(:,1) = S(:,1) + s;
S(:,1+k)=0;

% Use cardinality constraints to strengthen bound
n_fixed_one = nnz(p.lb(p.binary_variables)==1);
not_fixed = p.binary_variables(find(p.lb(p.binary_variables)==0 & p.ub(p.binary_variables)==1));
n_free = p.binarycardinality.up - n_fixed_one;
must_be_set = p.binarycardinality.up == p.binarycardinality.down;
cont_var = find(p.c);
c =  p.c(cont_var);
upper_bound = -inf;
% S(cont)*cont + rest >= 0
for i = 1:size(S,1)
    if S(i,1+cont_var) && ...
            isempty(setdiff(find(S(i,2:end)),[cont_var(:); p.binary_variables(:)]))
        % At most n_free of the remaining variables can be set to 1
        % Hence with q = sum of the k largest positive elements
        % S(1) + q + S(cont)*xc >= 0  
        % Our cost is f+c*xc and thus we can derive a lower bound on
        % the objective
        if sign(c) == sign(S(i,1+cont_var))          
            row = full(S(i,1+not_fixed));
            if ~must_be_set
                % At most n_free binaries are set.
                row = max(row,0);
            end
            row = sort(row,'descend');
            count = max(0,min(n_free,numel(row)));
            max_binary_contribution = sum(row(1:count));
            bound = p.f - c*(S(i,1) + max_binary_contribution)/S(i,1+cont_var);
            upper_bound = max(upper_bound,bound);
        end
    end
    if upper_bound > upper
        % We have detected that the bound is poor here
        break
    end
end
