function localatmost = postdetector(p)

newCuts = [];
localatmost.groups = {};
localatmost.bounds = [];
return
LPs = p.F_struc(1:p.K.f+p.K.l,:);
for j = 1:length(p.K.s)
    top = startofSDPCone(p.K);
    if p.K.s(j) >= 3
        Z = p.F_struc(top:top+p.K.s(j)^2-1,:);
        Zp = reshape(any(Z(:,2:end),2),p.K.s(j),[]);
        r = symamd(Zp);
        X = reshape(1:p.K.s(j)^2,p.K.s(j),[]);
        X = X(r,r);
        Z = Z(X(:),:);
        n = 3;
        X = spalloc(p.K.s(j),p.K.s(j),p.K.s(j)^2);
        X(1:n,1:n) = 1;       
        index0 = find(X);
        index = index0;
        corner = 0;
        for block = 1:p.K.s(j)-n
            dataBlock = Z(index,:);
            vars = find(any(dataBlock(:,2:end),1));
            L = p.lb(vars);
            U = p.ub(vars);
            combs = dec2decbin(0:2^length(vars)-1,length(vars))';
            combs = repmat(L,1,size(combs,2)) + repmat(U-L,1,size(combs,2)).*combs;
            combs = unique(combs','rows')';
            dataBlock = dataBlock(:,[1 1+vars]);
            working = ones(1,size(combs,2));
            for k = 1:size(combs,2)
                if min(eig(reshape(dataBlock*[1;combs(:,k)],3,3)))<-1e-5
                    working(k)=0;
                end
            end
            if any(working)
                remaining = combs(:,find(working));
                q = repmat(1e30,length(p.c),1);
                for k = 1:size(remaining,2)
                    qq = q;q(vars) = remaining(:,k);
                    m = LPs*[1;qq];
                    m = m(abs(m)<=1e10);
                    if length(m)>0
                    if any(m)<0
                        remaining(:,k);
                    end
                    end
                end                    
            end
           
            if length(vars)>=2 %&& max(sum(combs(:,find(working)),1))==1                
                localatmost.groups{end+1} = vars;
                localatmost.bounds(end+1) = max(sum(combs(:,find(working)),1));
            end                          
            index = index + p.K.s(j)+1;
        end
    end
end      