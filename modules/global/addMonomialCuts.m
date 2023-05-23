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
                    % Derive tangent cuts at the borders, and add
                    % a cut parallell to the upper cut. The point where
                    % this is located is given by M
                    % The cut y>=0 will be added elsewhere
                    df = (U^n-L^n)/(U-L);
                    M = sign(df)*(abs(df)/n)^(1/(n-1));
                    [Ax,Ay,b,K] = convexhullConvex(L,M,U,L^n,M^n,U^n,n*L^(n-1),n*M^(n-1),n*U^(n-1));                    
                    p_cut.F_struc(end+1:end+length(b),1) = b;
                    p_cut.F_struc(end-length(b)+1:end,1+monom_variable) = -Ax;
                    p_cut.F_struc(end-length(b)+1:end,1+monom_index) = -Ay;
                    p_cut.K.l = p_cut.K.l+length(b);
                else
                    % Odd case a bit trickier, as signs on bounds will 
                    % determine convexity 
                    if L >= 0
                        % Convex case. Add tangent cut as in even case
                        % We know df is non-negative here
                        df = (U^n-L^n)/(U-L);
                        M = (df/n)^(1/(n-1));
                        [Ax,Ay,b,K] = convexhullConvex(L,M,U,L^n,M^n,U^n,n*L^(n-1),n*M^(n-1),n*U^(n-1));                    
                        p_cut.F_struc(end+1:end+length(b),1) = b;
                        p_cut.F_struc(end-length(b)+1:end,1+monom_variable) = -Ax;
                        p_cut.F_struc(end-length(b)+1:end,1+monom_index) = -Ay;
                        p_cut.K.l = p_cut.K.l+length(b);
                    elseif U <= 0
                        % Concave case. Add tangent cut as in even case
                        % We know df is non-negative here and we are
                        % looking for the tangent in negative space
                        df = (U^n-L^n)/(U-L);
                        M = -(df/n)^(1/(n-1));
                        [Ax,Ay,b,K] = convexhullConcave(L,M,U,L^n,M^n,U^n,n*L^(n-1),n*M^(n-1),n*U^(n-1));                    
                        p_cut.F_struc(end+1:end+length(b),1) = b;
                        p_cut.F_struc(end-length(b)+1:end,1+monom_variable) = -Ax;
                        p_cut.F_struc(end-length(b)+1:end,1+monom_index) = -Ay;
                        p_cut.K.l = p_cut.K.l+length(b);
                    else
                        % Tricky case neither convex nor concave
                        % Line between lower bound and tangent intersection
                        r = zeros(1,n+1);r(1)=n-1;r(2)=-L*n;r(end)=L^n;
                        r = roots(r);                        
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