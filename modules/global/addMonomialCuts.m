function p = addMonomialCuts(p)

if any(p.originalModel.variabletype==3)
    singleMonomial =  sum(p.originalModel.monomtable |  p.originalModel.monomtable,2)==1;
    p_cut = emptyNumericalModel;
    monomials = find(p.originalModel.variabletype==3);
    for i = 1:length(monomials)
        monom_index = monomials(i);
        monom = p.originalModel.monomtable(monomials(i),:);
        monom_variable = find(monom);
        if length(monom_variable)==1
            n = monom(monom_variable);
            L = p.lb(monom_variable);
            U = p.ub(monom_variable);
            if ~isinf(L) && ~isinf(U)
                if even(n)
                    M = (L+U)/2;
                    if L <= 0 && U >= 0 && M~=0
                        if M < 0
                            [Ax,Ay,b,K] = convexhullConvex(L,M,0,U,L^n,M^n,0,U^n,n*L^(n-1),n*M^(n-1),0,n*U^(n-1));
                        else
                            [Ax,Ay,b,K] = convexhullConvex(L,0,M,U,L^n,0,M^n,U^n,n*L^(n-1),0,n*M^(n-1),n*U^(n-1));
                        end
                    else
                        [Ax,Ay,b,K] = convexhullConvex(L,M,U,L^n,M^n,U^n,n*L^(n-1),n*M^(n-1),n*U^(n-1));
                    end
                    p_cut.F_struc(end+1:end+length(b),1) = b;
                    p_cut.F_struc(end-length(b)+1:end,1+monom_variable) = -Ax;
                    p_cut.F_struc(end-length(b)+1:end,1+monom_index) = -Ay;
                    p_cut.K.l = p_cut.K.l+length(b);
                else
                    if p.lb(monom_variable)<0 &  p.ub(monom_variable)>0 & ~isinf(p.lb(monom_variable)) & ~isinf(p.ub(monom_variable))
                        
                        % Line between lower bound and tangent intersection
                        r = zeros(1,n+1);r(1)=n-1;r(2)=-L*n;r(end)=L^n;
                        r = roots(r);
                        %r = r(min(find(r==real(r))));
                        r = max(r(imag(r)==0));
                        if r >= U
                            r = U;
                            fprim = (U^n-L^n)/(U-L);
                        else
                            fprim = n*r^(n-1);
                        end
                        p_cut.F_struc(end+1,1) = -L^n+L*fprim;
                        p_cut.F_struc(end,1+monom_index)=1;
                        p_cut.F_struc(end,1+monom_variable)=-fprim;
                        p_cut.K.l = p_cut.K.l+1;
                        
                        % Line between upper bound and tangent intersection
                        r = zeros(1,n+1);r(1)=n-1;r(2)=-U*n;r(end)=U^n;
                        r = roots(r);
                        r = min(r(imag(r)==0));
                        if r <= L
                            r = L;
                            fprim = (U^n-L^n)/(U-L);
                        else
                            fprim = n*r^(n-1);
                        end
                        p_cut.F_struc(end+1,1) = U^n-U*fprim;
                        p_cut.F_struc(end,1+monom_index)=-1;
                        p_cut.F_struc(end,1+monom_variable)=fprim;
                        p_cut.K.l = p_cut.K.l+1;
                    end
                end
            end
        end
    end
    p = mergeNumericalModels(p,p_cut);
end