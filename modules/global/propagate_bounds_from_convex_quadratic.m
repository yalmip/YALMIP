function p = presolve_quadratic_psdbound(p)
% Look for r >= (or ==) (x - c)^T Q (x - c)
if p.K.f + p.K.l > 0    
    for i = 1:p.K.l + p.K.f
        rhs = p.F_struc(i,1);
        c   = -p.F_struc(i,2:end);
        if ~any(c(find(p.variabletype>2))) && all(c(find(p.variabletype==2)))
            [Q,c] = compileQuadratic(c,p,0);
            Q = Q(p.linears,p.linears);
            c = c(p.linears);
            [R,e] = chol(Q);
            if ~e && nnz(c)==0
                A = diag(diag(Q) - sum(abs(Q-diag(diag(Q))),2));
                [R,e] = chol(A);
                if ~e
                    U = rhs./diag(A).^.5;
                    L = -rhs./diag(A).^.5;
                    p.ub(p.linears) = min([p.ub(p.linears),U],[],2);
                    p.lb(p.linears) = max([p.lb(p.linears),L],[],2);
                end
            end
        end
    end
end
