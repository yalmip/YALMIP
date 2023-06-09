function [p_lp,infeasibility,infeasible_socp_cones] = add_socp_cut(p,p_lp,x,infeasibility)
infeasible_socp_cones = zeros(1,length(p.K.q));
% Only add these cuts if solver doesn't support SOCP cones
if ~p.solver.lower.constraint.inequalities.secondordercone.linear
    if any(p.K.q)
        % Add cuts
        top = startofSOCPCone(p.K);
        for i = 1:1:length(p.K.q)
            n = p.K.q(i);
            X = p.F_struc(top:top+n-1,:)*[1;x];
            X = [X(1) X(2:end)';X(2:end) eye(n-1)*X(1)];
            Y = randn(n,n);
            newcuts = 1;
            newF = zeros(n,size(p.F_struc,2));
            [d,v] = eig(X);
            infeasibility = min(infeasibility,min(diag(v)));
            dummy=[];
            newF = [];
            if infeasibility<0
                [ii,jj] = sort(diag(v));
                for m = jj(1:min(length(jj),p.options.cutsdp.cutlimit))'%find(diag(v<0))%1:1%length(v)
                    if v(m,m)<0
                        v1 = d(1,m);v2 = d(2:end,m);
                        newF = [newF;p.F_struc(top,:) + 2*v1*v2'*p.F_struc(top+1:top+n-1,:)];
                        newcuts = newcuts + 1;
                    end
                end
            end
            newF(abs(newF)<1e-12) = 0;
            keep= any(newF(:,2:end),2);
            newF = newF(keep,:);
            if size(newF,1)>0
                p_lp.F_struc = [p_lp.F_struc;newF];
                p_lp.K.l = p_lp.K.l + size(newF,1);
                [i,j] = sort(p_lp.F_struc*[1;x]);
            end
            top = top+n;
        end
    end
end