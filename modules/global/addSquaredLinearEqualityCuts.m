function p = addMultipliedEqualityCuts(p)

newRows = [];
for j = 1:p.K.f
    row = p.F_struc(j,:);
    used = find(row(2:end));
    if ~any(p.variabletype(used))
        % We only add these cuts on models where all terms are used
        % in objective (out of lazyness so we now all monomials exist...)
        if all(all(p.nonshiftedQP.Q(used,used)))
            newRow = spalloc(1,length(row),0);
            % Linear equality
            q = row(1 + used);
            Q = full(q'*q);Q = triu(Q + triu(Q)-diag(diag(Q)));
            newRow = spalloc(1,length(row),nnz(Q));
            newRow(1) = row(1)^2;
            for i = 1:length(used)
                for k = i:length(used)
                    xv = used(i);
                    yv = used(k);
                    zv = findrows(p.bilinears(:,2:end),[xv yv]);
                    if isempty(zv)
                        zv = findrows(p.bilinears(:,2:end),[yv xv]);
                        if isempty(zv)
                            %weird!
                            return
                        end  
                    end
                    zv = p.bilinears(zv,1);
                    newRow(1,zv + 1) = -Q(i,k);
                end
            end 
            newRows = [newRows;newRow];
        end        
    end   
end
if ~isempty(newRows)
    p.F_struc = [newRows;p.F_struc];
    p.K.f = p.K.f + size(newRows,1);    
end