function pout = propagate_bounds_from_arbitrary_quadratics(p)

pout = p;
if p.bilinears~=0
    F_struc = p.F_struc;
    
    p.F_struc = [-p.F_struc(1:p.K.f,:);p.F_struc];
    p.K.f=2*p.K.f;
    
    InequalityConstraintState = [p.EqualityConstraintState;p.EqualityConstraintState;p.InequalityConstraintState];
    
    if p.K.l+p.K.f>0
        quadratic_variables = find(p.bilinears(:,2) == p.bilinears(:,3));
        if ~isempty(quadratic_variables)
            quadratic_variables = p.bilinears(quadratic_variables,1);
            for i = 1:length(quadratic_variables)
                k = quadratic_variables(i);
                
                if (p.lb(k) < p.ub(k)-1e-4)                   
                    x = p.bilinears(p.bilinears(:,1)==k,2);% x^2
                    candidates = find((InequalityConstraintState==1) & p.F_struc(1:p.K.f+p.K.l,1+k))';
                    for j = candidates
                        a = p.F_struc(j,2:end);
                        b = p.F_struc(j,1);
                        aij = a(k);
                        
                        % rest + cij*x + (+-)x^2 >= 0
                        cij = a(x);
                        cij = cij/abs(aij);
                        b = b/abs(aij);
                        a = a/abs(aij);
                        %  aij = sign(aij);
                        
                        a(k) = 0;
                        a(x) = 0;
                        indNEG = find(a < 0);
                        indPOS = find(a > 0);
                        LB = p.lb;
                        UB = p.ub;
                        LB(k) = 0;
                        UB(k) = 0;
                        LB(x) = 0;
                        UB(x) = 0;
                        a(k) = 0;
                        a(x) = 0;
                        
                        if aij < 0
                            % rest + cij*x - x^2 >= 0
                            
                            % Derive upper bound on rest
                            rest = (b+a([indPOS(:);indNEG(:)])*[UB(indPOS);LB(indNEG)]);
                            % Write as rest >= (x - center)^2-(center)^2
                            center = cij/2;
                            radii2 = rest + (cij/2)^2;
                            if radii2 > 0                                
                                newUB = center + sqrt(radii2);
                                newLB = center - sqrt(radii2);
                                p.lb(x) = max(p.lb(x),newLB);
                                p.ub(x) = min(p.ub(x),newUB);
                            end
                        else
                            % rest + cij*x + x^2 >= 0
                            % (x+cij/2)^2 - (cij/2)^2 + rest >= 0
                            % (x+cij/2)^2 >= radii
                            % Lower bound on rest
                            rest = (b+a([indPOS(:);indNEG(:)])*[UB(indPOS);LB(indNEG)]);
                            center = -cij/2;
                            radii2 = -rest + (cij/2)^2;
                            if radii2 > 0
                                left  = center-sqrt(radii2);
                                right = center+sqrt(radii2);
                                if p.ub(x) < right
                                    p.ub(x) = min(p.ub(x),left);
                                end
                                if p.lb(x) > left
                                    p.lb(x) = max(p.lb(x),right);
                                end
                            end
                        end                                
                    end
                    
                elseif p.lb(k)<0
                    
                end
            end
        end
    end
    
    if p.K.l+p.K.f>0
        bilinear_variables = find(p.bilinears(:,2) ~= p.bilinears(:,3));
        if ~isempty(bilinear_variables)
            bilinear_variables = p.bilinears(bilinear_variables,1);
            for i = 1:length(bilinear_variables)
                k = bilinear_variables(i);
                if p.lb(k) >= -5000000000 & (p.lb(k) < p.ub(k)-1e-4)
                    x = p.bilinears(p.bilinears(:,1)==k,2);% x^2
                    y = p.bilinears(p.bilinears(:,1)==k,3);% x^2
                    candidates = find((InequalityConstraintState==1) & p.F_struc(1:p.K.f+p.K.l,1+k))';
                    for j = candidates
                        a = p.F_struc(j,2:end);
                        aij = a(k);
                        if aij > 0
                            indNEG = find(a < 0);
                            indPOS = find(a > 0);
                            LB = p.lb;
                            UB = p.ub;
                            LB(k) = 0;
                            UB(k) = 0;
                            a(k) = 0;
                            newLB = (-p.F_struc(j,1)-a([indPOS(:);indNEG(:)])*[UB(indPOS);LB(indNEG)])/aij;
                            p.lb(k) = max(p.lb(k),newLB);
                            
                        elseif aij < 0
                            indNEG = find(a < 0);
                            indPOS = find(a > 0);
                            LB = p.lb;
                            UB = p.ub;
                            LB(k) = 0;
                            UB(k) = 0;
                            a(k) = 0;
                            newUB = (p.F_struc(j,1)+a([indPOS(:);indNEG(:)])*[UB(indPOS);LB(indNEG)])/(-aij);
                            p.ub(k) = min(p.ub(k),newUB);
                        end
                    end
                    
                    if p.lb(k)>0 & p.lb(x)>0 & p.lb(y)>0
                        p.ub(x) = min(p.ub(x), p.ub(k)/p.lb(y));
                        p.ub(y) = min(p.ub(y), p.ub(k)/p.lb(x));
                        p.lb(x) = max(p.lb(x), p.lb(k)/p.ub(y));
                        p.lb(y) = max(p.lb(y), p.lb(k)/p.ub(x));
                    end
                    
                    
                elseif p.lb(k)<0
                    
                end
            end
        end
    end
end

if ~isequal([p.lb p.ub],[pout.lb pout.ub])
    quad_v = find(p.bilinears(:,2) == p.bilinears(:,3));
    quad_x = p.bilinears(quad_v,2);
    quad_v = p.bilinears(quad_v,1);
    if ~isempty(quad_v)
        % y = x^2, x>=0, y >= L means x >= sqrt(L)
        k = find(p.lb(quad_x)>=0 & p.lb(quad_v)>0);
        if ~isempty(k)           
            p.lb(quad_x(k)) = max(p.lb(quad_x(k)),sqrt(p.lb(quad_v(k))));
        end
    end
end

if ~isequal([p.lb p.ub],[pout.lb pout.ub])
    pout.changedbounds = 1;
end
pout.lb = p.lb;
pout.ub = p.ub;