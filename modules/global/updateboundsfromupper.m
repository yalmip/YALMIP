function  p = updateboundsfromupper(p,upper,ppoly);
if ~isinf(upper)
    LU = [p.lb p.ub];
    if nnz(p.c.*(p.ub-p.lb)) == 1 & nnz(p.Q)==0
        i = find(p.c.*(p.ub-p.lb));
        if p.c(i)>0
            p.ub(i) = min([p.ub(i) upper]);
        end
    end
    if ~isempty(p.bilinears) & nnz(p.Q)==0
        quad_v = find(p.bilinears(:,2) == p.bilinears(:,3));
        quad_x = p.bilinears(quad_v,2);
        quad_v = p.bilinears(quad_v,1);
        i = find(p.c(quad_v)>0);
        quad_v = quad_v(i);
        quad_x = quad_x(i);
        if ~isempty(quad_v) & ~any(isinf(p.lb(quad_x))) & ~any(isinf(p.ub(quad_x)))
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
                rhs = upper - p.f - (h(h<0)'*p.ub(h<0) + h(h>=0)'*p.lb(h>=0)) + xc'*D*xc;
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
    
end
