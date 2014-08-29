function  [oldlb,oldub] = presolve_probebinary(A,b,c,lb,ub,binary_variables);

goon = 1;

AL0A  = (A>0).*A;
AG0A  = (A<0).*A;

At = A';
%used = full(any(A(:,find(changed_bounds)),2));
isbinary  = ismembc(1:length(lb),binary_variables);
isinteger = ismembc(1:length(lb),binary_variables);

goon = all(lb<=ub);

oldub = ub;
oldlb = lb;

for i = 1:length(binary_variables)
        
    % Probe up
    lb = oldlb;
    ub = oldub;
    if lb(binary_variables(i))~=ub(binary_variables(i))
        lb(binary_variables(i)) = 1;   
        bi_dn = AL0A*lb+AG0A*ub;
        if any(bi_dn>b)
            oldlb(binary_variables(i)) = 0;   
            oldub(binary_variables(i)) = 0;   
        end
    end
    % Probe down
    lb = oldlb;
    ub = oldub;
    if lb(binary_variables(i))~=ub(binary_variables(i))
        ub(binary_variables(i)) = 0;   
        bi_dn = AL0A*lb+AG0A*ub;
        if any(bi_dn>b)
            oldlb(binary_variables(i)) = 1;   
            oldub(binary_variables(i)) = 1;   
        end
    end
end

fix_this = [];
fix_value = [];