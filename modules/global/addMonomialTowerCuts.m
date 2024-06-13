function p = addMonomialTowerCuts(p)

if any(p.originalModel.variabletype==3)
    singleMonomial =  sum(p.originalModel.monomtable |  p.originalModel.monomtable,2)==1;
    p_cut = emptyNumericalModel;
    monomials = find(p.originalModel.variabletype>=2 & p.originalModel.variabletype<=2);
    for i = 1:length(monomials)
        monom_index = monomials(i);
        monom = p.originalModel.monomtable(monomials(i),:);
        monom_variable = find(monom);
        if length(monom_variable)==1
            n = monom(monom_variable);
            L = p.lb(monom_variable);
            U = p.ub(monom_variable);
            if ~isinf(L) & ~isinf(U)
                if even(n)
                    candidates = find(singleMonomial & p.originalModel.monomtable(:,monom_variable)>n);
                    for k = candidates(:)'
                        m = p.originalModel.monomtable(k,monom_variable);
                        if even(m)
                            % We are looking at the graph of [x^n;x^m] and want to
                            % add a tangent cut i.e. lower bound by (dx^m/dx^n)
                            % f(x) >= f(x0) + (g(x) - g(x0))*df(x0)/dg(x0);
                            if L <= 0 && U >= 0
                                Lg = 0;Ug = max(L^n,U^n);
                            else
                                Lg = min(L^n,U^n);Ug = max(L^n,U^n);
                            end
                            x0 = ((Lg + Ug)/2)^(1/n);
                            y0 = x0^m;
                            df = (m/n)*x0^(m-n);
                            x1 = ((Lg))^(1/n);
                            x2 = ((Ug))^(1/n);
                            
                            [Ax,Ay,b,K] = convexhullConvex(Lg,((Lg + Ug)/2),Ug,Lg^(m/n),y0,Ug^(m/n),(m/n)*x1^(m-n),df,(m/n)*x2^(m-n));
                            p_cut.F_struc(end+1:end+length(b),1) = b;
                            p_cut.F_struc(end-length(b)+1:end,1+k) = -Ay;
                            p_cut.F_struc(end-length(b)+1:end,1+monom_index) = -Ax;
                            p_cut.K.l = p_cut.K.l+length(b);
                        end
                    end
                end
            end
        end
    end
    p = mergeNumericalModels(p,p_cut);
end