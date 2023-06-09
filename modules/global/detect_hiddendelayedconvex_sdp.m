function p = detect_hiddendelayedconvex_sdp(p)
p.hiddendelayedconvex.variable = [];
p.hiddendelayedconvex.modelup = {};
p.hiddendelayedconvex.modeldown = {};
if any(p.K.s)
    top = startofSDPCone(p.K);
    for i = 1:length(p.K.s)
        n = p.K.s(i);
        m = size(p.F_struc,2);
        digindex = top + find(speye(p.K.s(i)))-1;
        f = p.F_struc(digindex,:);
        possible_q = find((f(:,1)<0) & (sum(f|f,2)==2));
        possible_constant = find((f(:,1)>0) & (sum(f|f,2)==1));
        for j = possible_q(:)'
            q = find(f(j,2:end));
            if f(j,q+1) > 0 && p.variabletype(q)==2
                y = p.QuadraticsList(q);
                for k = possible_constant(:)'
                    % we have something like [y^2-w ?;? c]
                    w = -f(j,1);
                    c = f(k,1);
                    off_index = (j-1)*n + k;
                    r = p.F_struc(top + off_index-1,:);
                    vars = find(r(2:end));
                    if ~any(p.variabletype(vars))
                        % So [y^2-w r;r c]>=0
                        % y^2 >= w + r^2/c  dvs y >=norm([w^.5;r/sqrt(c)])
                        U = emptyNumericalModel;
                        U.F_struc = [zeros(1,m);
                            sqrt(w) zeros(1,m-1);
                            r/sqrt(c)];
                        U.F_struc(1,1+y)=1;
                        U.K.q = 3;
                        L = emptyNumericalModel;
                        L.F_struc = [zeros(1,m);
                            sqrt(w) zeros(1,m-1);
                            r/sqrt(c)];
                        L.F_struc(1,1+y)=-1;
                        L.K.q = 3;
                        p.hiddendelayedconvex.variable = [p.hiddendelayedconvex.variable [y;0;0]];
                        p.hiddendelayedconvex.modelup{end+1}=U;
                        p.hiddendelayedconvex.modeldown{end+1}=L;                            
                    end
                end
            end
        end
        top = top + n^2;
    end
end