function p = presolve_empty_rows(p,dontprint)

if p.feasible
    n = 0;
    if p.K.f > 0
        j = find(abs(p.F_struc(1:p.K.f,1))<eps & ~any(p.F_struc(1:p.K.f,2:end),2));
        p.F_struc(j,:)=[];
        p.K.f = p.K.f - length(j);
        n = n + length(j);
    end
    
    if p.K.l > 0
        j = find(abs(p.F_struc(startofLPCone(p.K):endOfLPCone(p.K),1)<eps) & ~any(p.F_struc(startofLPCone(p.K):endOfLPCone(p.K),2:end),2));
        p.F_struc(p.K.f + j,:)=[];
        p.K.l = p.K.l - length(j);
        n = n + length(j);
    end
    
    if p.options.verbose>1 && n > 0 & nargin < 2
        disp(['* Removed ' num2str(n) ' empty rows']);
    end
end