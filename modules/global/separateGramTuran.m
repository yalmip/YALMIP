function [cuts,maxViolation] = separateGramTuran(gram,x,nvars,objective)
% A saturated Gram block represents at most K distinct vector classes.
% For r vertices, balanced K-way classes minimize the number of equal pairs.
cuts = sparse(0,nvars+1);
maxViolation = 0;
for block = 1:numel(gram)
    g = gram{block};
    n = size(g.pair,1);
    d = g.diagonal;
    b = g.otherEntry;
    if ~g.saturated || ~g.exactPattern || d <= 0 || b >= 0 || ...
            nnz(triu(g.pair,1)) ~= n*(n-1)/2
        continue
    end
    pairVariables = unique(nonzeros(g.pair(:)));
    if nargin >= 4 && ...
            nnz(objective(pairVariables)) < 0.5*numel(pairVariables)
        continue
    end
    K = floor(1-d/b);
    if K < 2 || K >= n || d+K*b >= -1e-8*(1+abs(d))
        continue
    end
    weight = zeros(n);
    weight(g.pair > 0) = x(g.pair(g.pair > 0));
    if g.equalValue == 0
        weight = 1-weight;
    end
    sizes = unique([min(2*K,n),n]);
    blockCuts = sparse(0,nvars+1);
    blockViolations = [];
    for r = sizes
        q = floor(r/K);
        s = mod(r,K);
        rhs = s*q*(q+1)/2 + (K-s)*q*(q-1)/2;
        if rhs == 0
            continue
        end
        candidates = zeros(0,r);
        violations = [];
        seeds = 1:n;
        if r == n
            seeds = 1;
        end
        for seed = seeds
            if r == n
                selected = 1:n;
            else
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
            end
            total = sum(sum(triu(weight(selected,selected),1)));
            violation = rhs-total;
            maxViolation = max(maxViolation,violation);
            if violation > 1e-6
                candidates(end+1,:) = selected; %#ok<AGROW>
                violations(end+1,1) = violation; %#ok<AGROW>
            end
        end
        if isempty(candidates)
            continue
        end
        [candidates,keep] = unique(candidates,'rows');
        violations = violations(keep);
        for k = 1:size(candidates,1)
            subset = candidates(k,:);
            edgeVars = nonzeros(triu(g.pair(subset,subset),1));
            counts = ones(numel(edgeVars),1);
            if g.equalValue == 1
                row = sparse(counts,1+edgeVars,counts,1,nvars+1);
                row(1) = -rhs;
            else
                row = sparse(counts,1+edgeVars,-counts,1,nvars+1);
                row(1) = numel(edgeVars)-rhs;
            end
            blockCuts = [blockCuts;row]; %#ok<AGROW>
            blockViolations(end+1,1) = violations(k); %#ok<AGROW>
        end
    end
    if ~isempty(blockViolations)
        [~,order] = sort(blockViolations,'descend');
        cuts = [cuts;blockCuts(order(1:min(20,numel(order))),:)]; %#ok<AGROW>
    end
end
end
