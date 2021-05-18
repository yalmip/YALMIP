function res = computeconstraintresiduals(p,x)

res= [];
if ~isempty(p.F_struc)
    vecres = p.F_struc*[1;x];
    if any(p.K.f)
        res = -abs(vecres(1:p.K.f));
    end
    if any(p.K.l)
        res = [res;vecres(p.K.f+1:p.K.f+p.K.l)];
    end    
    if any(p.K.q)
        top = startofSOCPCone(p.K)
        for i = 1:length(p.K.q)
            n = p.K.q(i);
            X = vecres(top:top+n-1);top = top+n;
            res = [res;X(1)-norm(full(X(2:end)))];
        end
    end
    if any(p.K.s)
        top = startofSDPCone(p.K);
        for i = 1:length(p.K.s)
            n = p.K.s(i);
            X = reshape(vecres(top:top+n^2-1),n,n);top = top+n^2;
            res = [res;min(eig(X))];
        end
    end
end
