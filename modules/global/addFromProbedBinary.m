function p = addFromProbedBinary(p)

% For a nonlinear operator, the envelope model is tighter on a smaller
% domain. Search for (binary d implies x >= L) and then compute bounds on 
% f(x) on that domain, and add implies(d, L <= f(x) <= U)
if ~isempty(p.binary_variables) && ~isempty(p.F_struc)
    newCuts = zeros(0,size(p.F_struc,2));
    b = p.F_struc(:,1);
    A = p.F_struc(:,2:end);
    D = A(:,p.binary_variables);
    A(:,p.binary_variables) = 0;
    probedBoundsL = [];
    probedBoundsU = [];
    for i = p.K.f + 1:p.K.f + p.K.l
        [~,ix,valx] = find(A(i,:));
        if length(ix) == 1
            [~,id,vald] = find(D(i,:));
            if length(id) == 1
                if valx > 0
                    % b + valx*x + vald*d >= 0
                    % valx*x >= (-b-vald*x)
                    % d implies x >= (-b-vald)/valx
                    activated_lower = (-b(i)-vald)/valx;
                    probedBoundsL(end+1,1) = id;
                    probedBoundsL(end,2) = ix;
                    probedBoundsL(end,3) = activated_lower;
                elseif valx < 0
                    % b - (-valx)*x + vald*d >= 0
                    % b + vald*d >= -valx*x
                    % d implies (b+vald/(-valx)
                    activated_upper = (b(i)+vald)/(-valx);
                    probedBoundsU(end+1,1) = id;
                    probedBoundsU(end,2) = ix;
                    probedBoundsU(end,3) = activated_upper;
                end
            end
        end
    end
    if ~isempty(probedBoundsL) || ~isempty(probedBoundsU)
        for i = 1:length(p.binary_variables)
            iL = find(probedBoundsL(:,1) == i);
            iU = find(probedBoundsU(:,1) == i);
            if ~isempty(iL) || ~isempty(iU)
                ptemp = p;
                if ~isempty(iL)
                    for k = iL
                        ptemp.lb(probedBoundsL(k,2)) = probedBoundsL(k,3);
                    end
                end
                if ~isempty(iU)
                    for k = iU
                        ptemp.ub(probedBoundsU(k,2)) = probedBoundsU(k,3);
                    end
                end
                ptemp.lb(p.binary_variables(i)) = 1;
                pprobed = update_monomial_bounds(ptemp);
                pprobed = propagate_bounds_from_evaluations(pprobed);              
                sL = find(pprobed.lb > ptemp.lb+1e-6);
                sU = find(pprobed.ub < ptemp.ub-1e-6);
                for j = 1:length(sL)
                    % {d implies x(s) >= q}
                    % x(s) >= q*d + L*(1-d)
                    % -L + x(s) + d*(L-q) >= 0
                    if ~ismember( sL(j),p.binary_variables)   % {d implies x(s) <= q}
                        q = pprobed.lb(sL(j));
                        L = p.lb(sL(j));
                        newCuts(end+1,1) = -L;
                        newCuts(end,1+sL(j)) = 1;
                        newCuts(end,1+p.binary_variables(i)) = L-q;
                    end
                end
                for j = 1:length(sU)
                    if ~ismember( sU(j),p.binary_variables)   % {d implies x(s) <= q}
                        % x(s) <= q*d + U*(1-d)
                        % U - x(s) + d*(q-U) >= 0
                        q = pprobed.ub(sU(j));
                        U = p.ub(sU(j));
                        newCuts(end+1,1) = U;
                        newCuts(end,1+sU(j)) = -1;
                        newCuts(end,1+p.binary_variables(i)) = q-U;
                    end
                end
            end
        end
    end
    if nnz(newCuts)>0
        p.F_struc = [p.F_struc(1:p.K.f+p.K.l,:);
            newCuts;
            p.F_struc(p.K.f+p.K.l+1:end,:)];
        p.K.l = p.K.l + size(newCuts,1);
    end
end