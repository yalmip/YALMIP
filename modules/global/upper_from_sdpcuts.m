function upper_bound = upper_from_sdpcuts(p,sdpCuts,upper)
% In models here we have a simple unit continuous objective
% we can use sdp cuts collected on the way to quickly derive
% a crude bound
% If this bound is larger than upper bound, there is no
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
    if S(i,1+cont_var)
        % At most n_free of the remaining variables can be set to 1
        % Hence with q = sum of the k largest positive elements
        % S(1) + q + S(cont)*xc >= 0  
        % Our cost is c*xc and thus we can derive an upper bound on
        % the objective
        if sign(c) == sign(S(i,1+cont_var))          
            if n_free == 0
                upper_bound = max(upper_bound,-c*(S(i,1))./S(i,1+cont_var));
            else
                row = S(i,1+p.binary_variables);
                row = S(i,1+not_fixed);
                if ~must_be_set
                    % We only set variables to 1 when positive
                    row = max(row,0);              
                end              
                try
                    if isempty(row)
                        upper_bound = max(upper_bound,-c*(S(i,1))./S(i,1+cont_var));              
                    else
                        upper_bound = max(upper_bound,-c*(S(i,1) + sumk(row,n_free))./S(i,1+cont_var));              
                    end
                catch
                    1
                end
            end
        end
    end
    if upper_bound > upper
        % We have detected that the bound is poor here
        break
    end
end