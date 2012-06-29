function feasible = checkfeasiblefast(p,x,tol)

feasible = 0;
if ~(all(x < p.ub+tol) & all(x > p.lb-tol));
    return
end

if ~isempty(p.F_struc)
    vecres = p.F_struc*[1;x];

    if p.K.f>0
        if any( -abs(vecres(1:p.K.f)) < - tol)
            return
        end
    end

    if p.K.l>0
        if any(vecres(p.K.f+1:p.K.f+p.K.l) <-tol)
            return
        end
    end

    if p.K.q(1)>0
        top = 1+p.K.f+p.K.l;
        for i = 1:length(p.K.q)
            n = p.K.q(i);
            X = vecres(top:top+n-1);top = top+n;
            if any(X(1)-norm(full(X(2:end))) < -tol)
                return
            end
        end
    end

    if p.K.s(1)>0
        top = 1+p.K.f+p.K.l+p.K.q;
        for i = 1:length(p.K.s)
            n = p.K.s(i);
            X = reshape(vecres(top:top+n^2-1),n,n);top = top+n^2;
            X = (X+X')/2;
            if min(real(eig(X))) < -tol
                return
            end
        end
    end
end
feasible = 1;