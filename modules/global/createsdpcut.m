function p = createsdpcut(p,x)
if any(p.K.s) 
    % Cuts are added from violated relaxations
    % Might improve bound propagation.
    top = startofSDPCone(p.K);
    newcuts = 1;
    newF = [];    
    newSOCP = [];
    for i = 1:length(p.K.s)
        n = p.K.s(i);
        X = p.F_struc(top:top+n^2-1,:)*[1;x];
        X = full(reshape(X,n,n));
        [d,v] = eig(X);
        if any(diag(v)<0)
             for m = 1:length(v)
                if v(m,m)<0
                    for j = 1:length(x)+1
                        newF(newcuts,j)= d(:,m)'*reshape(p.F_struc(top:top+n^2-1,j),n,n)*d(:,m);
                    end
                    newF(newcuts,1)=newF(newcuts,1)+1e-6;
                    newcuts = newcuts + 1;
                    if size(p.lpcuts,1)>0
                        dist = p.lpcuts*newF(newcuts-1,:)'/(newF(newcuts-1,:)*newF(newcuts-1,:)');
                        if any(abs(dist-1)<1e-3)
                            newF = newF(1:end-1,:);
                            newcuts = newcuts - 1;
                        end
                    end
                end
            end
        end        
        top = top+n^2;
    end

    if ~isempty(newF)
        m = size(newF,2);        
        p.lpcuts = [newF;p.lpcuts];
        p.cutState = [ones(size(newF,1),1);p.cutState];          
    end
    if ~isempty(newSOCP)
        p.socpcuts.F_struc = [p.socpcuts.F_struc;newSOCP];
        p.socpcuts.K.q = [p.socpcuts.K.q 3];    
    end
end