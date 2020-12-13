function [res,X,detailed] = constraint_residuals(p,x)
res= [];
X = [];
detailed.f = [];
detailed.l = [];
detailed.q = [];
detailed.s = [];
detailed.b = [];
if p.K.f>0    
    detailed.f = -abs(p.F_struc(1:p.K.f,:)*[1;x]);
    res = detailed.f;
end
if p.K.l>0
    detailed.l = p.F_struc(p.K.f+1:p.K.f+p.K.l,:)*[1;x];
    res = [res;detailed.l];
end
if nnz(p.K.q)>0
    z = p.F_struc(1 + p.K.f+p.K.l:p.K.f+p.K.l+sum(p.K.q),:)*[1;x];
    top = 1;
    for i = 1:length(p.K.q)
        zi = z(top:top + p.K.q(i)-1);top = top + p.K.q(i);
        detailed.q = [detailed.q;zi(1)^2-zi(2:end)'*zi(2:end)];
    end
    res = [res;detailed.q];
end

if (length(p.K.s)>1) | p.K.s>0
    top = 1+p.K.f+p.K.l + sum(p.K.q);
    for i = 1:length(p.K.s)
        n = p.K.s(i);
        X{i} = p.F_struc(top:top+n^2-1,:)*[1;x];top = top+n^2;
        X{i} = reshape(X{i},n,n);
        detailed.s = [detailed.s;min(eig(X{i}))];
    end
    res = [res;detailed.s];
end
detailed.b = min([p.ub-x;x-p.lb]);
res = [res;detailed.b];