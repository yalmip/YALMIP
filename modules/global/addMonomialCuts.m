function p = addMonomialCuts(p)

if any(p.originalModel.variabletype==3)
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
                    
            if even(n)
%                 M = (L+U)/2;
%                 [Ax,Ay,b,K] = convexhullConvex(L,M,U,L^n,M^n,U^n,n*L^(n-1),n*M^(n-1),n*U^(n-1));
%                 p.F_struc(end+1:end+length(b),1) = b;
%                 p.F_struc(end-length(b)+1:end,1+monom_variable) = -Ax;
%                 p.F_struc(end-length(b)+1:end,1+monom_index) = -Ay;
%                 p.K.l = p.K.l+length(b);
            else
                if p.lb(monom_variable)<0 &  p.ub(monom_variable)>0
                   
                    % Tangent at x = lower bound
                    if L<=0                        
                        p_cut.F_struc(end+1,1) = L^n-n*L^n;
                        p_cut.F_struc(end,1+monom_index)=-1;
                        p_cut.F_struc(end,1+monom_variable)=n*L^(n-1);
                        p_cut.K.l = p_cut.K.l+1;
                    end
                     
                    % Tangent at x = upper bound
                    if U >= 0
                        p_cut.F_struc(end+1,1) = -(U^n-n*U^n);
                        p_cut.F_struc(end,1+monom_index)= 1;
                        p_cut.F_struc(end,1+monom_variable)=-n*U^(n-1);
                        p_cut.K.l =  p_cut.K.l+1;
                    end
                    
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
                    
%                     % Line between upper bound and tangent intersection
                    r = zeros(1,n+1);r(1)=n-1;r(2)=-U*n;r(end)=U^n;
                    r = roots(r);
                    %r = r(min(find(r==real(r))));
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
    p = mergeNumericalModels(p,p_cut);
end