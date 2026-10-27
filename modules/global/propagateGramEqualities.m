function [p,changed] = propagateGramEqualities(p,gram)
% Equality in a PSD Gram block identifies two rows, so incident pair
% variables have equal values. This only applies when one pair value
% exactly equals the constant diagonal.
changed = false;
if ~p.feasible
    return
end
for block = 1:numel(gram)
    if ~gram{block}.saturated
        continue
    end
    pair = gram{block}.pair;
    equalValue = gram{block}.equalValue;
    n = size(pair,1);
    parent = 1:n;
    for j = 2:n
        for i = 1:j-1
            v = pair(i,j);
            if v && p.lb(v) >= equalValue-1e-10 && ...
                    p.ub(v) <= equalValue+1e-10
                parent(findRoot(parent,j)) = findRoot(parent,i);
            end
        end
    end
    roots = zeros(1,n);
    for i = 1:n
        roots(i) = findRoot(parent,i);
    end
    components = unique(roots);
    if numel(components) == n
        continue
    end
    members = cell(1,numel(components));
    for a = 1:numel(components)
        members{a} = find(roots == components(a));
    end
    for a = 1:numel(components)
        membersA = members{a};
        if numel(membersA) == 1
            continue
        end
        internal = unique(nonzeros(pair(membersA,membersA)));
        if ~isempty(internal)
            oldLower = p.lb(internal);
            oldUpper = p.ub(internal);
            p.lb(internal) = max(oldLower,equalValue);
            p.ub(internal) = min(oldUpper,equalValue);
            changed = changed || any(p.lb(internal) > oldLower) || ...
                any(p.ub(internal) < oldUpper);
            if any(p.lb(internal) > p.ub(internal)+1e-12)
                p.feasible = 0;
                return
            end
        end
        for b = 1:numel(components)
            if a == b
                continue
            end
            membersB = members{b};
            edges = unique(nonzeros(pair(membersA,membersB)));
            if numel(edges) < 2
                continue
            end
            lower = max(p.lb(edges));
            upper = min(p.ub(edges));
            if lower > upper+1e-12
                p.feasible = 0;
                return
            end
            changed = changed || any(p.lb(edges) < lower) || ...
                any(p.ub(edges) > upper);
            p.lb(edges) = max(p.lb(edges),lower);
            p.ub(edges) = min(p.ub(edges),upper);
        end
    end
end
end

function root = findRoot(parent,i)
root = i;
while parent(root) ~= root
    root = parent(root);
end
end
