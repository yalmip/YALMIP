function p = createsdpcut(p,x)
if p.K.s > 0
    top = p.K.f+p.K.l+sum(p.K.q)+1;
    newcuts = 1;
    newF = [];
    for i = 1:length(p.K.s)
        n = p.K.s(i);
        X = p.F_struc(top:top+n^2-1,:)*[1;x];
        X = full(reshape(X,n,n));
        [d,v] = eig(X);
        for m = 1:length(v)
            if v(m,m)<0
                for j = 1:length(x)+1;
                    newF(newcuts,j)= d(:,m)'*reshape(p.F_struc(top:top+n^2-1,j),n,n)*d(:,m);
                end
                % max(abs(newF(:,2:end)),[],2)
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
        top = top+n^2;
    end

    if ~isempty(newF)
        % Don't keep all
        m = size(newF,2);
        %  size(p.lpcuts)
        p.lpcuts = [newF;p.lpcuts];
        p.cutState = [ones(size(newF,1),1);p.cutState];
     %   violations = p.lpcuts*[1;x];
     %   p.lpcuts = p.lpcuts(violations<0.1,:);
     %   p.cutState = p.cutState(violations<0.1);
        
        if size(p.lpcuts,1)>15*m
            %disp('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!');
            %violations = p.lpcuts*[1;x];
            %[i,j] = sort(violations);
            %p.lpcuts = p.lpcuts(j(1:15*m),:);
            %p.cutState = lpcuts = p.lpcuts(j(1:15*m),:);
            %p.lpcuts = p.lpcuts(end-15*m+1:end,:);
        end
    end
end