function p = detect3x3SymmetryGroups(p)
   
if ~isempty(p.binary_variables)
    groups = {};
    for j = 1:length(p.K.s)
        if p.K.s(j) >= 3 && p.semidefinite{j}.purebinary
            n = 3;
            X = spalloc(p.K.s(j),p.K.s(j),p.K.s(j)^2);
            X(1:n,1:n) = 1;
            index0 = find(X);
            index = index0;
            corner = 0;
            for block = 1:p.K.s(j)-n+1
                dataBlock = p.semidefinite{j}.F_struc(index,:);
                used = find(any(dataBlock,1));
                dataBlock = dataBlock(:,used);
                if used(1) > 1
                    dataBlock = [zeros(size(dataBlock,1),1) dataBlock];
                end                
                v = used;v = v(v>1)-1;
                if all(p.isbinary(v))
                    if isempty(groups)
                        groups{1}.dataBlock = dataBlock;
                        groups{1}.variables{1} = used;
                    else
                        found = 0;
                        for k = 1:length(groups)
                            if isequal(length(used),length(groups{1}.variables{1}))
                                if isequal(groups{1}.dataBlock,dataBlock)
                                    found = 1;
                                    groups{1}.variables{end+1} = used;
                                else
                                    % TODO: Look for simple scaled versions
                                    % s = groups{1}.dataBlock(:)\dataBlock(:);
                                    % norm(s*groups{1}.dataBlock-dataBlock,inf)
                                end
                            end
                        end
                        if ~found
                            groups{end+1}.dataBlock = dataBlock;
                            groups{end}.variables{1} = used;
                        end
                    end
                end
                index = index + p.K.s(j)+1;
            end
        end
    end
    for i = 1:length(groups)
        if 1%length(groups{i}.variables) > 1
            keep(i) = 1;
        else
            keep(i) = 0;
        end
    end
    if length(groups) > 0
        groups = {groups{find(keep)}};
        if length(groups) > 0
            for i = 1:length(groups)
                for j = 1:length(groups{i}.variables)
                    v = groups{i}.variables{j};
                    v = v(v>1)-1;
                    groups{i}.variables{j} = v;
                end
            end
        end
    end
else
    groups = {};
end
p.sdpsymmetry = groups;
