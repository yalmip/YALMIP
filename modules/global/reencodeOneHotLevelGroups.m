function [p,groups] = reencodeOneHotLevelGroups(p,candidates)
%REENCODEONEHOTLEVELGROUPS Branch on cumulative choices in one-hot groups.
% For a one-hot binary group x_1,...,x_k, use threshold binaries
% y_j=sum_{i>=j}x_i. The inverse x_i=y_i-y_{i+1} is applied to every
% original model coefficient, and y_j>=y_{j+1} enforces the same domain.
% Unlike coefficient aggregation, this coordinate change is exact even
% when nominally proportional SDP coefficients differ by file rounding.
groups = {};
if numel(candidates)<4 || isempty(p.K.s) || p.K.f || ...
        nnz(p.Q) || any(p.variabletype) || ~isempty(p.evalMap) || ...
        ~isempty(p.nonlinear)
    return
end
candidate = intersect(candidates(:),p.binary_variables(:));
candidate = candidate(p.lb(candidate)==0 & p.ub(candidate)==1);
if numel(candidate)<4
    return
end
linear = p.F_struc(1:p.K.l,:);
coneRows = p.K.l+1:size(p.F_struc,1);
used = false(numel(p.c),1);
for row = 1:p.K.l
    if abs(full(linear(row,1))-1)>1e-12
        continue
    end
    members = find(linear(row,2:end));
    k = numel(members);
    if k<2 || any(used(members)) || ~all(ismember(members,candidate)) || ...
            any(abs(full(linear(row,1+members))+1)>1e-12)
        continue
    end
    % Choose increasing structural levels to make threshold branches
    % split the ordered SDP effects. This is a performance gate only;
    % all subsequent coefficients use the exact original values.
    A = p.F_struc(coneRows,1+members);
    scales = full(max(abs(A),[],1));
    if any(scales<=1e-9)
        continue
    end
    [scales,order] = sort(scales);
    members = members(order);
    levels = 1:k;
    if max(abs(scales/scales(1)-levels))>1e-6 || ...
            max(abs(p.c(members(:))-levels(:)*p.c(members(1))))> ...
            1e-7*max(1,max(abs(p.c(members))))
        continue
    end
    base = full(A(:,order(1)));
    ordered = full(A(:,order));
    if max(max(abs(ordered-base*levels)))> ...
            1e-7*max(1,max(abs(ordered(:))))
        continue
    end
    groups{end+1} = members(:)'; %#ok<AGROW>
    used(members) = true;
end
if isempty(groups)
    return
end
added = sum(cellfun(@numel,groups))-numel(groups);
newRows = sparse(added,numel(p.c)+1);
next = 0;
for groupIndex = 1:numel(groups)
    members = groups{groupIndex};
    k = numel(members);
    originalColumns = p.F_struc(:,1+members);
    originalCost = p.c(members);
    p.F_struc(:,1+members(2:end)) = ...
        originalColumns(:,2:end)-originalColumns(:,1:end-1);
    p.c(members(2:end)) = originalCost(2:end)-originalCost(1:end-1);
    if ~isempty(p.x0)
        originalGuess = p.x0(members);
        p.x0(members) = flipud(cumsum(flipud(originalGuess(:))));
    end
    for level = 1:k-1
        next = next+1;
        newRows(next,1+members(level)) = 1;
        newRows(next,1+members(level+1)) = -1;
    end
end
oldLinearCount = p.K.l;
p.F_struc = [p.F_struc(1:oldLinearCount,:);newRows; ...
    p.F_struc(oldLinearCount+1:end,:)];
p.K.l = oldLinearCount+added;
p.rowtype = [p.rowtype(:)' zeros(1,added)];
end
