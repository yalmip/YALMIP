function p = presolve_empty_rows(p)

if p.K.f > 0
    j = find(abs(p.F_struc(1:p.K.f,1))<eps & ~any(p.F_struc(1:p.K.f,2:end),2));
    p.F_struc(j,:)=[];
    p.K.f = p.K.f - length(j);
end

if p.K.l > 0
    j = find(abs(p.F_struc(startofLPCone(p.K):endOfLPCone(p.K),1)<eps) & ~any(p.F_struc(startofLPCone(p.K):endOfLPCone(p.K),2:end),2));
    p.F_struc(p.K.f + j,:)=[];
    p.K.l = p.K.l - length(j);
end