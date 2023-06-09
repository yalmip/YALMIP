function p = presolve_dualreductions(p)

if p.feasible && nnz(p.Q)==0
    n = 0;
    enters_linearly = (sum(p.monomtable | p.monomtable,1) == 1) & (sum(p.monomtable,1) == 1);
    for i = 1:length(p.c)
        if enters_linearly(i)
            [row,~,val] = find(p.F_struc(:,1+i));
            if ~isempty(row) && all(row > p.K.f) && all(row <= p.K.f + p.K.l)
                % This variable is only involved in LP constraints
                s = sign(val);
                if all(s(1) == s)
                    % and enters with same sign
                    % If we don't care about variable, just
                    if p.c(i) == 0
                        if s(1)>0 && ~isinf(p.ub(i))
                            p.lb(i) = p.ub(i);
                            n = n+1;
                        elseif ~isinf(p.lb(i))
                            p.ub(i) = p.lb(i);
                            n = n+1;
                        end
                    elseif p.c(i)>0 && s(1)<0
                        % We want to make it small, and from feasibility
                        % that is optimal too
                        if ~isinf(p.lb(i))
                            p.ub(i) = p.lb(i);
                            n = n+1;
                        end
                    elseif  p.c(i)<0 && s(1)>0
                        % We want to make it large, and from feasibility
                        % that is optimal too.
                        if ~isinf(p.ub(i))
                            p.lb(i) = p.ub(i);
                            n = n+1;
                        end
                    end
                end
            end
        end
    end
    if n > 0 && p.options.verbose > 1
        disp(['* ' num2str(n) ' variables fixed by dual reduction']);
    end
end