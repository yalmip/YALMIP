function  p = propagate_bounds_from_upper(p,upper)
if nargin == 1
    upper = p.upper;
end
if upper ~= fix(upper)
    upper = upper + 1e-7;
end
if ~isinf(upper)
    LU = [p.lb p.ub];
    % Simple objective f + c_i*x(i)
    if nnz(p.c) == 1 & nnz(p.Q)==0
        i = find(p.c);
        if p.c(i) > 0
            % We are minimizing x(i), since an upper bound is UPPER
            % this means c(i)*x(i) has to be < UPPER in optimal solution
            p.ub(i) = min([p.ub(i) (upper-p.f)/p.c(i)]);
        elseif p.c(i) < 0
            % We are maximizing x(i), since an lower bound is -UPPER
            % this means x(i) has to be > -UPPER in optimal solution
            p.lb(i) = max([p.lb(i) -(upper-p.f)/abs(p.c(i))]);
        end            
    end
    if nnz(p.Q)==0
        % Very basic, simply propagate from f + sum c_i m_i(x) <= upper
        c = p.c;
        i = find(c);
        for j = i(:)'
            cc = c;cc(j)=0;
            ii = find(cc);
           % ii = setdiff(i,j);% Much slower!
            if ~isempty(ii)
                if p.c(j) > 0
                    obound = upper - p.f;
                    ii_pos = ii(find(p.c(ii) > 0));
                    ii_neg = ii(find(p.c(ii) < 0));
                    obound = obound-sum(p.c(ii_pos).*p.lb(ii_pos));
                    obound = obound-sum(p.c(ii_neg).*p.ub(ii_neg));
                    p.ub(j) = min(p.ub(j),obound/p.c(j));
                else
                    obound = p.f - upper;
                    ii_pos = ii(find(p.c(ii) > 0));
                    ii_neg = ii(find(p.c(ii) < 0));
                    obound = obound+sum(p.c(ii_pos).*p.lb(ii_pos));
                    obound = obound+sum(p.c(ii_neg).*p.ub(ii_neg));
                    p.lb(j) = max(p.lb(j),obound/-p.c(j));                   
                end
            end
        end
    end
    if ~isempty(p.bilinears) & nnz(p.Q)==0
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
        i = find(p.c(quad_v)>0);
        quad_v = quad_v(i);
        quad_x = quad_x(i);
        if ~isempty(quad_v) %& ~any(isinf(p.lb(quad_x))) & ~any(isinf(p.ub(quad_x)))
            % x'*D*x + c'*x + f + h(x) < U
            D = diag(p.c(quad_v));
            if all(diag(D)>0)
                c = p.c(quad_x);
                h = p.c;
                h(quad_x) = 0;
                h(quad_v) = 0;
                % Complete
                % (x + D^-0.5*c/2)'*D*(x + D^-0.5*c/2) - c'*D*c/4 + f + h(x) < U
                xc = -(D^-1)*c/2;
                %rhs = upper - p.f - (h(h<0)'*p.ub(h<0) + h(h>=0)'*p.lb(h>=0)) + xc'*D*xc;
                rhs = upper - p.f + xc'*D*xc;
                ix = find(h<0);
                if ~isempty(ix)
                    rhs = rhs - h(ix)'*p.ub(ix);
                end
                ix = find(h>0);
                if ~isempty(ix)
                    rhs = rhs - h(ix)'*p.lb(ix);
                end
                if rhs>0
                    D = D/rhs;
                    % (x-xc)'*D*(x-xc) <= 1
                    for i = 1:length(quad_x)
                        p.lb(quad_x(i)) = max(p.lb(quad_x(i)),xc(i) - 1/sqrt(D(i,i)));
                        p.ub(quad_x(i)) = min(p.ub(quad_x(i)),xc(i) + 1/sqrt(D(i,i)));
                    end
                end
            end
        end
    end
    if ~isequal(LU,[p.lb p.ub])
        p.changedbounds = 1;
    end
    
    % Some initial code for using inverse objective to derive bounds.
    if nnz(p.Q)==0 && nnz(p.c)==1
        [pos,~,val] = find(p.c);
        if val == -1            
            fi = find(p.evalVariables == pos);
            if ~isempty(fi)
                % We're maximizing f(x)
                if ~isempty(p.evalMap{fi}.properties.inverse)
                    if strcmp(p.evalMap{fi}.properties.definiteness,'positive')
                        if strcmp(p.evalMap{fi}.properties.monotonicity,'increasing')
                            if length(p.evalMap{fi}.arg)==2                               
                                lower = -upper;
                                p.lb(p.evalMap{fi}.variableIndex) = max([p.lb(p.evalMap{fi}.variableIndex) p.evalMap{fi}.properties.inverse(lower)]);
                            end
                        end
                    end
                end
            end
        end
        if val == 1
            fi = find(p.evalVariables == pos);
            if ~isempty(fi)
                % We're minimizing f(x)
                if ~isempty(p.evalMap{fi}.properties.inverse)
                    if strcmp(p.evalMap{fi}.properties.definiteness,'positive')
                        if strcmp(p.evalMap{fi}.properties.monotonicity,'increasing')
                            if length(p.evalMap{fi}.arg)==2                                                              
                                p.ub(p.evalMap{fi}.variableIndex) = min([p.ub(p.evalMap{fi}.variableIndex) p.evalMap{fi}.properties.inverse(upper)]);
                            end
                        end
                    end
                end
            end
        end
    end
    % Numerical issues easily propagates, so widen weird close to feasible box
    p = widenSmashedBox(p);
    if any(p.lb > p.ub + 1e-7)
        p.feasible = 0;
    end
end