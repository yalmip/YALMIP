function p = addMonomialCuts(p)

if any(p.originalModel.variabletype==3)
    monomials = find(p.originalModel.variabletype==3);
    for i = 1:length(monomials)
        monom_index = monomials(i);
        monom = p.originalModel.monomtable(monomials(i),:);
        monom_variable = find(monom);
        if length(monom_variable)==1
            n = monom(monom_variable);
            if ~even(n)
                if p.lb(monom_variable)<0 &  p.ub(monom_variable)>0
                    L = p.lb(monom_variable);
                    U = p.ub(monom_variable);
                    
                    % Tangent at x = lower bound
                    if L<=0
                        p.F_struc(end+1,1) = L^n-n*L^n;
                        p.F_struc(end,1+monom_index)=-1;
                        p.F_struc(end,1+monom_variable)=n*L^(n-1);
                        p.K.l = p.K.l+1;
                    end
                     
                    % Tangent at x = upper bound
                    if U >= 0
                        p.F_struc(end+1,1) = -(U^n-n*U^n);
                        p.F_struc(end,1+monom_index)= 1;
                        p.F_struc(end,1+monom_variable)=-n*U^(n-1);
                        p.K.l = p.K.l+1;
                    end
                    
                    % Line between lower bound and tangent intersection
                    r = zeros(1,n+1);r(1)=n-1;r(2)=-L*n;r(end)=L^n;
                    r = roots(r);
                    r = r(min(find(r==real(r))));
                    if r >= U
                        r = U;
                        fprim = (U^n-L^n)/(U-L);
                    else
                        fprim = n*r^(n-1);
                    end
                    p.F_struc(end+1,1) = -L^n+L*fprim;
                    p.F_struc(end,1+monom_index)=1;
                    p.F_struc(end,1+monom_variable)=-fprim;
                    p.K.l = p.K.l+1;
                    
                    % Line between upper bound and tangent intersection
                    r = zeros(1,n+1);r(1)=n-1;r(2)=-U*n;r(end)=U^n;
                    r = roots(r);
                    r = r(min(find(r==real(r))));
                    if r <= L
                        r = L;
                        fprim = (U^n-L^n)/(U-L);
                    else
                        fprim = n*r^(n-1);
                    end
                    p.F_struc(end+1,1) = U^n-U*fprim;
                    p.F_struc(end,1+monom_index)=-1;
                    p.F_struc(end,1+monom_variable)=fprim;
                    p.K.l = p.K.l+1;
                end
            end
        end
    end
end