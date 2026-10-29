function [cuts,maxViolation] = separateGramCliques(gram,x,nvars,objective)
% Greedy separation of pigeonhole cuts certified by Gram minors.
% The detector has certified that equal-value edges are transitive. If
% d+(r-1)*b<0, an r-vertex all-other-value pattern is not PSD, so at
% least one pair must have the equal value. Exact saturation is unnecessary.
cuts = sparse(0,nvars+1);
maxViolation = 0;
for block = 1:numel(gram)
    g = gram{block};
    n = size(g.pair,1);
    d = g.diagonal;
    b = g.otherEntry;
    if ~g.exactPattern || d <= 0 || b >= 0 || ...
            nnz(triu(g.pair,1)) ~= n*(n-1)/2
        continue
    end
    if nargin >= 4
        pairVariables = unique(nonzeros(g.pair(:)));
        if nnz(objective(pairVariables)) < 0.5*numel(pairVariables)
            continue
        end
    end
    % SDPA decimals can put -d/b just below an integer. Find the first
    % cardinality with a strictly negative all-other eigenvalue instead
    % of flooring a numerically rounded ratio.
    margin = 1e-5*(1+abs(d)+abs(b));
    K = 2;
    while K < n && d+K*b >= -margin
        K = K+1;
    end
    r = K+1;
    if r > n || d+K*b >= -margin
        continue
    end
    weight = zeros(n);
    weight(g.pair > 0) = x(g.pair(g.pair > 0));
    if g.equalValue == 0
        weight = 1-weight;
    end
    candidate = zeros(0,r);
    violations = [];
    for seed = 1:n
        selected = seed;
        available = true(1,n);
        available(seed) = false;
        while numel(selected) < r
            choices = find(available);
            scores = sum(weight(selected,choices),1);
            [~,best] = min(scores);
            vertex = choices(best);
            selected(end+1) = vertex; %#ok<AGROW>
            available(vertex) = false;
        end
        for pass = 1:2
            improved = false;
            for position = 1:r
                old = selected(position);
                others = selected([1:position-1 position+1:r]);
                choices = find(available);
                gains = sum(weight(others,choices),1) - ...
                    sum(weight(others,old));
                [gain,best] = min(gains);
                if ~isempty(gain) && gain < -1e-9
                    selected(position) = choices(best);
                    available(old) = true;
                    available(choices(best)) = false;
                    improved = true;
                end
            end
            if ~improved
                break
            end
        end
        selected = sort(selected);
        total = sum(sum(triu(weight(selected,selected),1)));
        violation = 1-total;
        maxViolation = max(maxViolation,violation);
        if violation > 1e-6
            candidate(end+1,:) = selected; %#ok<AGROW>
            violations(end+1,1) = violation; %#ok<AGROW>
        end
    end
    if isempty(candidate)
        continue
    end
    [candidate,keep] = unique(candidate,'rows');
    violations = violations(keep);
    [~,order] = sort(violations,'descend');
    candidate = candidate(order(1:min(20,numel(order))),:);
    blockCuts = sparse(size(candidate,1),nvars+1);
    for k = 1:size(candidate,1)
        subset = candidate(k,:);
        edgeVars = nonzeros(triu(g.pair(subset,subset),1));
        counts = ones(numel(edgeVars),1);
        if g.equalValue == 1
            blockCuts(k,:) = sparse(counts,1+edgeVars,counts,1,nvars+1);
            blockCuts(k,1) = -1;
        else
            blockCuts(k,:) = sparse(counts,1+edgeVars,-counts,1,nvars+1);
            blockCuts(k,1) = numel(edgeVars)-1;
        end
    end
    cuts = [cuts;blockCuts]; %#ok<AGROW>
end
end
