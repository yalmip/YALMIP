function p = detect_quadratic_disjoints_generalized(p)

p.quadraticdisjoints = [];
if any(p.K.l)
    top = startofLPCone(p.K);
    for i = 1:p.K.l
        row = p.F_struc(top,:);
        if row(1)<0
            [Q, c] =  compileQuadratic(row(2:end),p,0); 
            Q = -Q;
            c = -c;
            if ~any(c)
                used = find(any(Q,2));
                Q = Q(used,used);
                if ~(min(eig(Q))>=0)
                    for r = 1:length(used)
                        s = -row(1);
                        k = -Q(r,r);
                        a = Q(:,r);
                        a(r,:) = [];
                        Qsub=Q;Qsub(r,:)=[];Qsub(:,r)=[];
                        if min(eig(Qsub))>0
                            w = sqrt(s/(k+a'*inv(Qsub)*a));
                            p.quadraticdisjoints = [p.quadraticdisjoints [used(r);(w);-(w)]];
                        end
                    end
                end
            end            
        end
        top = top + 1;
    end
end
