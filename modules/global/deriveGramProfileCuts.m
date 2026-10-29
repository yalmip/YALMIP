function [cuts,profiles] = deriveGramProfileCuts(p,gram)
% Derive class-size cuts from a certified complete binary Gram block.
% The Gram data must come from initGramTriangleSeparator. Each equal-value
% pair belongs to one equivalence class, and a negative all-different minor
% bounds the number of classes. Exact degree rows then bound each class size.
% Each recognized linear row bounds the number of equal neighbors of one
% vertex. The resulting class sizes are enumerated by a small integer DP.
nvars = size(p.F_struc,2)-1;
cuts = sparse(0,nvars+1);
profiles = struct('n',{},'K',{},'rawDegree',{},'degree',{}, ...
    'equalPairs',{},'classCounts',{},'addedRows',{},'reason',{});
if ~isfield(p.K,'f') || ~isfield(p.K,'l') || p.K.l <= 0
    return
end
linear = p.F_struc(p.K.f+(1:p.K.l),:);
for block = 1:numel(gram)
    g = gram{block};
    n = size(g.pair,1);
    record = struct('n',n,'K',NaN,'rawDegree',[],'degree',[], ...
        'equalPairs',[],'classCounts',[],'addedRows',0,'reason','');
    edges = nonzeros(triu(g.pair,1));
    nedge = n*(n-1)/2;
    if ~g.exactPattern || numel(edges) ~= nedge || ...
            numel(unique(edges)) ~= nedge || ...
            g.diagonal <= 0 || g.otherEntry >= 0 || ...
            ~all(p.isbinary(edges))
        record.reason = 'Gram certificate or pair map unavailable';
        profiles(end+1) = record; %#ok<AGROW>
        continue
    end
    margin = 1e-5*(1+abs(g.diagonal)+abs(g.otherEntry));
    K = 2;
    while K < n && g.diagonal+K*g.otherEntry >= -margin
        K = K+1;
    end
    if K >= n || g.diagonal+K*g.otherEntry >= -margin
        record.reason = 'No strict class-count witness';
        profiles(end+1) = record; %#ok<AGROW>
        continue
    end
    if K > 12
        record.reason = 'Class-count DP cap';
        profiles(end+1) = record; %#ok<AGROW>
        continue
    end
    record.K = K;

    % For each edge variable store its two vertex endpoints. A row is used
    % only if its exact support is one complete incident star, its pair
    % coefficients are uniformly +1 or -1, and its constant is an integer.
    [vi,vj] = find(triu(g.pair,1));
    edgeids = full(g.pair(sub2ind([n n],vi,vj)));
    incidence = sparse([vi;vj],[edgeids;edgeids],1,n,nvars);
    lower = -inf(n,1);
    upper = inf(n,1);
    for rowIndex = 1:size(linear,1)
        row = linear(rowIndex,:);
        [~,columns,coefficients] = find(row(:,2:end));
        if numel(columns) ~= n-1 || ...
                ~all(coefficients == coefficients(1)) || ...
                abs(coefficients(1)) ~= 1
            continue
        end
        constant = full(row(1));
        if ~isfinite(constant) || constant ~= round(constant)
            continue
        end
        incident = find(sum(incidence(:,columns),2) == n-1);
        if numel(incident) ~= 1
            continue
        end
        vertex = incident(1);
        if g.equalValue == 1
            if coefficients(1) == 1
                lower(vertex) = max(lower(vertex),-constant);
            else
                upper(vertex) = min(upper(vertex),constant);
            end
        else
            if coefficients(1) == 1
                upper(vertex) = min(upper(vertex),n-1+constant);
            else
                lower(vertex) = max(lower(vertex),n-1-constant);
            end
        end
    end
    if any(~isfinite(lower)) || any(~isfinite(upper)) || ...
            any(lower ~= lower(1)) || any(upper ~= upper(1)) || ...
            lower(1) < 0 || upper(1) > n-1 || lower(1) > upper(1)
        record.reason = 'No uniform finite degree interval';
        profiles(end+1) = record; %#ok<AGROW>
        continue
    end
    L = lower(1);
    U = upper(1);
    record.rawDegree = [L U];
    sizes = L+1:U+1;
    reachable = false(n+1,K+1);
    minimum = inf(n+1,K+1);
    maximum = -inf(n+1,K+1);
    reachable(1,1) = true;
    minimum(1,1) = 0;
    maximum(1,1) = 0;
    for classes = 0:K-1
        for used = 0:n
            if ~reachable(used+1,classes+1)
                continue
            end
            for sizeValue = sizes
                target = used+sizeValue;
                if target > n
                    continue
                end
                equalPairs = sizeValue*(sizeValue-1)/2;
                reachable(target+1,classes+2) = true;
                minimum(target+1,classes+2) = min( ...
                    minimum(target+1,classes+2), ...
                    minimum(used+1,classes+1)+equalPairs);
                maximum(target+1,classes+2) = max( ...
                    maximum(target+1,classes+2), ...
                    maximum(used+1,classes+1)+equalPairs);
            end
        end
    end
    classCounts = find(reachable(n+1,2:end));
    if isempty(classCounts)
        record.reason = 'No feasible class-size profile';
        profiles(end+1) = record; %#ok<AGROW>
        continue
    end
    minPairs = min(minimum(n+1,classCounts+1));
    maxPairs = max(maximum(n+1,classCounts+1));
    possibleSizes = [];
    for sizeValue = sizes
        if any(reachable(n-sizeValue+1,1:K))
            possibleSizes(end+1) = sizeValue; %#ok<AGROW>
        end
    end
    degreeRange = [min(possibleSizes)-1 max(possibleSizes)-1];
    record.degree = degreeRange;
    record.equalPairs = [minPairs maxPairs];
    record.classCounts = classCounts;
    newCuts = sparse(0,nvars+1);
    % Summing the existing per-vertex degree rows gives
    % n*L/2 <= E <= n*U/2. Do not emit the weaker integral bounds.
    usefulLower = minPairs > ceil(n*L/2);
    usefulUpper = maxPairs < floor(n*U/2);
    if g.equalValue == 1
        if usefulLower
            newCuts(end+1,:) = sparse(ones(1,nedge+1),[1 1+edges'], ...
                [-minPairs ones(1,nedge)],1,nvars+1); %#ok<AGROW>
        end
        if usefulUpper
            newCuts(end+1,:) = sparse(ones(1,nedge+1),[1 1+edges'], ...
                [maxPairs -ones(1,nedge)],1,nvars+1); %#ok<AGROW>
        end
    else
        if usefulLower
            newCuts(end+1,:) = sparse(ones(1,nedge+1),[1 1+edges'], ...
                [nedge-minPairs -ones(1,nedge)],1,nvars+1); %#ok<AGROW>
        end
        if usefulUpper
            newCuts(end+1,:) = sparse(ones(1,nedge+1),[1 1+edges'], ...
                [maxPairs-nedge ones(1,nedge)],1,nvars+1); %#ok<AGROW>
        end
    end
    for vertex = 1:n
        incidentEdges = nonzeros(g.pair(vertex,:));
        if degreeRange(1) > L
            if g.equalValue == 1
                constant = -degreeRange(1);
                signValue = 1;
            else
                constant = n-1-degreeRange(1);
                signValue = -1;
            end
            newCuts(end+1,:) = sparse(ones(1,n),[1 1+incidentEdges'], ...
                [constant signValue*ones(1,n-1)],1,nvars+1); %#ok<AGROW>
        end
        if degreeRange(2) < U
            if g.equalValue == 1
                constant = degreeRange(2);
                signValue = -1;
            else
                constant = degreeRange(2)-(n-1);
                signValue = 1;
            end
            newCuts(end+1,:) = sparse(ones(1,n),[1 1+incidentEdges'], ...
                [constant signValue*ones(1,n-1)],1,nvars+1); %#ok<AGROW>
        end
    end
    record.addedRows = size(newCuts,1);
    record.reason = 'ok';
    cuts = [cuts;newCuts]; %#ok<AGROW>
    profiles(end+1) = record; %#ok<AGROW>
end
end
