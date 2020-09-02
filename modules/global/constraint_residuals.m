function [res,X] = constraint_residuals(p,x)
res= [];
X = [];
if p.K.f>0
    res = -abs(p.F_struc(1:p.K.f,:)*[1;x]);
end
if p.K.l>0
    res = [res;p.F_struc(p.K.f+1:p.K.f+p.K.l,:)*[1;x]];
end
if nnz(p.K.q)>0
    z = p.F_struc(1 + p.K.f+p.K.l:p.K.f+p.K.l+sum(p.K.q),:)*[1;x];
    top = 1;
    for i = 1:length(p.K.q)
        zi = z(top:top + p.K.q(i)-1);top = top + p.K.q(i);
        res = [res;zi(1)^2-zi(2:end)'*zi(2:end)];
    end
end

if (length(p.K.s)>1) | p.K.s>0
    top = 1+p.K.f+p.K.l + sum(p.K.q);
    for i = 1:length(p.K.s)
        n = p.K.s(i);
        X{i} = p.F_struc(top:top+n^2-1,:)*[1;x];top = top+n^2;
        X{i} = reshape(X{i},n,n);
        try
            res = [res;min(eig(X{i}))];
        catch
            1
        end
    end
end
res = [res;min([p.ub-x;x-p.lb])];