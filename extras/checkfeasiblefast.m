function feasible = checkfeasiblefast(p,x,tol)

feasible = 0;
if ~(all(x < p.ub+tol) && all(x > p.lb-tol))
    return
end

if ~isempty(p.F_struc)
    vecres = p.F_struc*[1;x];

    if any(p.K.f)
        if any(-abs(vecres(1:p.K.f)) < - tol)
            return
        end
    end

    if any(p.K.l)
        if any(vecres(p.K.f+1:p.K.f+p.K.l) <-tol)
            return
        end
    end

    if any(p.K.q)
        top = startofSOCPCone(p.K);
        for i = 1:length(p.K.q)
            n = p.K.q(i);
            X = vecres(top:top+n-1);top = top+n;
            if X(1)-norm(full(X(2:end))) < -tol
                return
            end
        end
    end
    
    if any(p.K.e)
        top = startofEXPCone(p.K);
        for i = 1:p.K.e           ;
            X = vecres(top:top+2);top = top+3;
            if X(2)==0 && X(3) < -tol
                return
            elseif X(3) - X(2)*exp(X(1)/X(2)) < -tol
                return
            end
        end
    end
    
    if any(p.K.p)
        top = startofPOWCone(p.K);
        for i = 1:length(p.K.p)
            X = vecres(top:top+p.K.p(i)-1);
            a = X(end);
            top = top+p.K.p(i);
            if X(1)^a*X(2)^(1-a) - norm(X(3:end-1)) < -tol
                return
            end
        end
    end

    if any(p.K.s)
        top = startofSDPCone(p.K);
        for i = 1:length(p.K.s)
            n = p.K.s(i);
            X = reshape(vecres(top:top+n^2-1),n,n);top = top+n^2;
            X = (X+X')/2 + tol*eye(n);
            [~,fail] = chol(X);
            if fail
                return
            end
        end
    end
end
feasible = 1;