function vars = decide_branch_variables(p)

if size(p.bilinears,1)==0
    if any(p.K.s)
        if any(p.K.s>p.K.rank)
            vars = p.linears;
            return
        end
    end
end

if p.options.bmibnb.lowrank==0
    nonlinear = find(~(sum(p.monomtable~=0,2)==1 & sum(p.monomtable,2)==1));
    vars      =  find(sum(abs(full(p.monomtable(nonlinear,:))),1));
    if ~isempty(p.evalVariables)
        temp = [];
        for i = 1:length(p.evalMap)
            temp = [temp p.evalMap{i}.variableIndex];
        end
        vars = union(vars,temp);
    end
else
    % Try to find a bi-partite structure
    pool1 = p.bilinears(1,2);
    pool2 = p.bilinears(1,3);
    
    for i = 2:size(p.bilinears,1)
        v1 = p.bilinears(i,2);
        v2 = p.bilinears(i,3);
        if v1==v2
            % We are fucked
            pool1 = [pool1 v1];
            pool2 = [pool2 v2];
        else
            if ismember(v1,pool1)
                pool2 = [pool2 v2];
            elseif ismember(v1,pool2)
                pool1 = [pool1 v2];
            elseif ismember(v2,pool1)
                pool2 = [pool2 v1];
            elseif ismember(v2,pool2)
                pool1 = [pool1 v1];
            else
                % No member yet
                pool1 = [pool1 v1];
                pool2 = [pool2 v2];
            end
        end
    end
    pool1 = unique(pool1);
    pool2 = unique(pool2);
    if isempty(intersect(pool1,pool2))
        if length(pool1)<=length(pool2)
            vars = pool1;
        else
            vars = pool2;
        end
    else
        nonlinear = find(~(sum(p.monomtable~=0,2)==1 & sum(p.monomtable,2)==1));
        vars =  find(sum(abs(full(p.monomtable(nonlinear,:))),1));
    end
end
