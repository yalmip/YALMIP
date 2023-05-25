function p = detect_cliques(p)

p.isclique = zeros(1,p.K.f+p.K.l);
for i = 1:p.K.f+p.K.l
    [~,idx,vals] = find(p.F_struc(i,:));
    if vals(1) == 1 && length(idx)>1 && idx(1) == 1
        idx = idx(2:end)-1;
        if all(p.isbinary(idx))
            if all(vals(2:end)==-1)
                p.isclique(i)=1;
            end
        end
    end
end
