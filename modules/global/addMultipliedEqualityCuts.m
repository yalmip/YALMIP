function p = addMultipliedEqualityCuts(p)

newRows = [];
for i =  p.linears
    for j = 1:p.K.f
        row = p.F_struc(j,:);
        used = find(row(2:end));
        if ~any(p.variabletype(used))
            % Linear constraint
            ok = 1;
            for k = 1:length(used)
                find1 = findrows(p.bilinears(:,2:end),[i used(k)]);
                find2 = findrows(p.bilinears(:,2:end),[used(k) i]);
                if (isempty(find1) & isempty(find2))
                    ok = 0;
                    break
                else
                    monom(k) = unique([find1 find2]);
                end
            end
            if ok
                monom = p.bilinears(monom,1);
                row(1+monom) = row(1+used);
                row(1+used) = 0;
                row(1+i) = row(1);
                row(1)=0;
                newRows = [newRows;row];
            end
        end                    
    end
end
if ~isempty(newRows)
    p.F_struc = [newRows;p.F_struc];
    p.K.f = p.K.f + size(newRows,1);
end