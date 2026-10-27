function bound = sdpCutLagrangianBound(p,sdpCuts,maxCuts,target)
%SDPCUTLAGRANGIANBOUND Lower bound from valid SDP eigenvector cuts.
% The caller supplies a model with a linear objective supported on binary
% variables. For each cut b + a'*z >= 0 and multiplier lambda >= 0,
% c'*z >= -lambda*b + (c-lambda*a)'*z. Minimize the right-hand side over
% this node's binary bounds. The maximum
% of these lower bounds is safe for pruning before the SDP node solve.

binary = p.binary_variables(:);
fixed = p.lb(binary) == p.ub(binary);
free = ~fixed;

bound = -inf;
c = full(p.c(binary));
firstCut = max(1,sdpCuts.K.l-maxCuts+1);
for cut = firstCut:sdpCuts.K.l
    row = sdpCuts.F_struc(cut,:);
    if ~all(ismember(find(row(2:end)),binary))
        continue
    end
    a = full(row(1+binary)');
    ratios = c(free)./a(free);
    ratios = sort(unique(ratios(isfinite(ratios) & ratios > 0)));
    if numel(ratios) > 5
        locations = unique(round(linspace(1,numel(ratios),5)));
        ratios = ratios(locations);
    end
    % The cut coefficients come from a numerical eigenvector. Relax its
    % constant by more than the coefficient cleanup in add_one_sdp_cut.
    cutSlack = 1e-8*(1+full(sum(abs(row))));
    for lambda = [0;ratios(:)]'
        reducedCost = c-lambda*a;
        constant = p.f-lambda*(row(1)+cutSlack) + ...
            reducedCost(fixed)'*p.lb(binary(fixed));
        candidate = constant + sum(min(0,reducedCost(free)));
        if isfinite(candidate)
            bound = max(bound,candidate);
        end
        if bound > target
            return
        end
    end
end
