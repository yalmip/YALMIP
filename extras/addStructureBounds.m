function [F_struc,K] = addStructureBounds(F_struc,K,ub,lb)
%ADDBOUNDS Internal function to add variable bounds to an F_struc

% Upper bounds
finite_bounds = find(~isinf(ub));
if ~isempty(finite_bounds)
    p = length(finite_bounds);
    n = length(ub);
    if ~isempty(F_struc)
        F_struc = [F_struc(1:1:K.f,:);[ub(finite_bounds) -sparse(1:p,finite_bounds,ones(p,1),p,n)];F_struc(K.f+1:end,:)];
    else
        F_struc = [ub(finite_bounds) -sparse(1:p,finite_bounds,ones(p,1),p,n)];
    end
    K.l = K.l+p;
end

% Lower bounds
finite_bounds = find(~isinf(lb));
if ~isempty(finite_bounds)   
    p = length(finite_bounds);
    n = length(ub);
    if ~isempty(F_struc)
        F_struc = [F_struc(1:1:K.f,:);[-lb(finite_bounds) sparse(1:p,finite_bounds,ones(p,1),p,n)];F_struc(K.f+1:end,:)];
    else
        F_struc = [[-lb(finite_bounds) sparse(1:p,finite_bounds,ones(p,1),p,n)];F_struc(K.f+1:end,:)];   
    end
    K.l = K.l+p;
end



