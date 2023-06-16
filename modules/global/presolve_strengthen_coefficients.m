function p = presolve_strengthen_coefficients(p)

if p.feasible % && isempty(p.nonlinear)
    if p.K.f>0
        Aeq = -p.F_struc(1:p.K.f,2:end);
        beq = p.F_struc(1:p.K.f,1);
        A = [Aeq;-Aeq];
        b = [beq;-beq];
        [p.lb,p.ub,redundant,pss] = tightenbounds(A,b,p.lb,p.ub,p.integer_variables,p.binary_variables,ones(length(p.lb),1));
    end
    pss=[];
    if p.K.l>0        
        A = -p.F_struc(1+p.K.f:p.K.f+p.K.l,2:end);
        b = p.F_struc(1+p.K.f:p.K.f+p.K.l,1);        
        if any(isinf(p.lb) | isinf(p.ub))
            % FIXME Added afterwards, should be integrated
            p = presolve_infs(A,b,p);
        end        
        [p.lb,p.ub,redundant,pss] = tightenbounds(A,b,p.lb,p.ub,p.integer_variables,p.binary_variables,ones(length(p.lb),1));
        if ~isempty(redundant)
            pss.AL0A(redundant,:)=[];
            pss.AG0A(redundant,:)=[];
            p.F_struc(p.K.f+redundant,:)=[];
            p.K.l = p.K.l - length(redundant);
        end
    end
end

function p = presolve_infs(A,b,p)
problematic = isinf(p.lb) | isinf(p.ub);
% Sort rows to contain few bad initially
% FIXME: Update problematic w.r.t direction
Abad = A(:,problematic);
[~,r] = sort(sum(Abad | Abad,2),'ascend');
r = r(r>0);
A = A(r,:);
b = b(r);
for i = 1:size(A,1)    
    row = A(i,:);
    [~,loc,val] = find(row);
    % If there is only one unbounded on this row, we can bound it
    if nnz(problematic(loc)) == 1
        bad = loc(find(problematic(loc)));
        L = p.lb;
        U = p.ub;
        L(bad) = 0;
        U(bad) = 0;
        rest_min = (row.*(row>0)*L + row.*(row<0)*U);
        rest_max = (row.*(row>0)*U + row.*(row<0)*L);
        if row(bad)>0
            p.ub(bad) = min(p.ub(bad),(b(i)-rest_min)/row(bad));
        else
            p.lb(bad) = max(p.lb(bad),(rest_min-b(i))/(-row(bad)));
        end
        problematic = isinf(p.lb) | isinf(p.ub);
    end
end
    

