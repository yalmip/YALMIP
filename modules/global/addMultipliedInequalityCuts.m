function p = addMultipliedInequalityCuts(p)

newRows = [];
for i =  p.linears
    if p.lb(i) >= -1000 | p.ub(i) <= 1000 % Avoid numerically tricky stuff
        for j = 1+p.K.f:p.K.f+p.K.l
            row = p.F_struc(j,:);
            used = find(row(2:end));
            if ~any(p.variabletype(used))
                % Linear constraint
                ok = 1;
                monom = zeros(1,length(used));
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
                    row0 = row;
                    if p.lb(i) >= -1000
                        row(1+monom) = row(1+used);
                        row(1+used) = 0;
                        row(1+i) = row(1);
                        row(1)=0;
                        newRows = [newRows;row-p.lb(i)*row0];
                    end
                    if p.ub(i) <= 1000
                        row = row0;
                        row(1+monom) = row(1+used);
                        row(1+used) = 0;
                        row(1+i) = row(1);
                        row(1)=0;
                        newRows = [newRows;-row+p.ub(i)*row0];
                    end
                end
            end
        end
    end
end
if ~isempty(newRows)
    p.F_struc = [p.F_struc;newRows];
    p.K.l = p.K.l + size(newRows,1);
end