function res = constraint_residuals(p,x)
res= [];
if p.K.f>0
    res = -abs(p.F_struc(1:p.K.f,:)*[1;x]);
end
if p.K.l>0
    res = [res;p.F_struc(p.K.f+1:p.K.f+p.K.l,:)*[1;x]];
end
if (length(p.K.s)>1) | p.K.s>0
    top = 1+p.K.f+p.K.l;
    for i = 1:length(p.K.s)
        n = p.K.s(i);
        X = p.F_struc(top:top+n^2-1,:)*[1;x];top = top+n^2;
        X = reshape(X,n,n);
        try
        res = [res;min(eig(X))];
        catch
            1
        end
    end
end
res = [res;min([p.ub-x;x-p.lb])];