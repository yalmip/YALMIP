function [prune,bound,info,candidate] = bnbDualCompressedScreen(p,upper)
%BNBDUALCOMPRESSEDSCREEN Bound and branch using a smaller auxiliary SDP.
% Its dual supplies a PSD direction; the bound uses the original node data.
% Its decoded primal is used only to choose a branch, never as an incumbent.
prune = false;
bound = -inf;
candidate = [];
info = struct('reason','inapplicable','compressedTime',0,'lpTime',0);
if ~isfinite(upper) || isempty(p.oneHotLevelGroups) || p.K.f || ...
        numel(p.K.s)~=1 || any(p.K.q) || nnz(p.Q) || ...
        any(p.variabletype) || ~isempty(p.evalMap) || ...
        ~isempty(p.semicont_variables) || ...
        ~strcmp(p.solver.lower.call,'callmosek')
    return
end
n = numel(p.c);
groups = p.oneHotLevelGroups;
coneRows = p.K.l+1:size(p.F_struc,1);
if numel(coneRows) ~= p.K.s(1)^2
    return
end
q = p;
q.options.saveduals = 1;
q.options.savesolverinput = 0;
q.options.savesolveroutput = 0;
q.x0 = [];
q.solver = p.solver.lower;
q.integer_variables = [];
q.binary_variables = [];
q.semicont_variables = [];
groupMask = false(n,1);
keepLinear = true(p.K.l,1);
for g = 1:numel(groups)
    members = groups{g}(:);
    if numel(members)<2 || any(groupMask(members)) || ...
            any(p.lb(members)<0) || any(p.ub(members)>1) || ...
            any(p.lb(members)~=round(p.lb(members))) || ...
            any(p.ub(members)~=round(p.ub(members))) || ...
            any(p.lb(members)>p.ub(members))
        info.reason = 'invalid group bounds';
        return
    end
    groupMask(members) = true;
    first = members(1);
    C = p.F_struc(coneRows,1+members);
    if nnz(C)==0
        info.reason = 'group absent from SDP';
        return
    end
    scale = max(1,max(abs(nonzeros(C))));
    if max(max(abs(C-C(:,1))))>1e-7*scale || ...
            max(abs(p.c(members)-p.c(first)))> ...
            1e-7*max(1,max(abs(p.c(members))))
        info.reason = 'weakly collinear group';
        return
    end
    lowerLevels = find(p.lb(members)>1-1e-10);
    upperLevels = find(p.ub(members)<1e-10);
    if isempty(lowerLevels), low = 0; else, low = max(lowerLevels); end
    if isempty(upperLevels)
        high = numel(members);
    else
        high = min(upperLevels)-1;
    end
    if low>high
        info.reason = 'empty group interval';
        return
    end
    q.lb(first) = low;
    q.ub(first) = high;
    q.lb(members(2:end)) = 0;
    q.ub(members(2:end)) = 0;
    q.c(members(2:end)) = 0;
    q.F_struc(coneRows,1+members(2:end)) = 0;
    linearColumns = p.F_struc(1:p.K.l,1+members);
    differences = linearColumns-linearColumns(:,1);
    linearScale = max(1,full(max(abs(linearColumns),[],2)));
    keepLinear = keepLinear & ...
        full(max(abs(differences),[],2))<=1e-9*linearScale;
    q.F_struc(1:p.K.l,1+members(2:end)) = 0;
end
% Retain rows whose group coefficients agree across cumulative thresholds.
% Others are dropped. This temporary model need not be a relaxation: its
% dual is used only as a PSD test matrix.
q.F_struc = [q.F_struc(find(keepLinear),:);q.F_struc(coneRows,:)];
q.K.l = nnz(keepLinear);
q.rowtype = zeros(1,q.K.l);
compressed = feval(p.solver.lower.call,q);
info.compressedTime = compressed.solvertime;
if compressed.problem~=0 || ~isfield(compressed,'Dual') || ...
        numel(compressed.Dual) < q.K.l+numel(coneRows)
    info.reason = 'compressed solve has no SDP dual';
    return
end
dimension = p.K.s(1);
% call_mosek_dual appends variable-bound rows before the SDP dual block.
% The SDP block is therefore at the end of output.Dual.
Z = reshape(full(compressed.Dual(end-numel(coneRows)+1:end)), ...
    dimension,dimension);
if any(~isfinite(Z(:)))
    info.reason = 'nonfinite SDP dual';
    return
end
Z = (Z+Z')/2;
[V,D] = eig(Z);
B = diag(sqrt(max(diag(D),0)))*V';
Z = B'*B;
lp = p;
lp.F_struc = p.F_struc(1:p.K.l,:);
lp.K.s = [];
lp.c = p.c-p.F_struc(coneRows,2:end)'*Z(:);
lp.f = p.f-p.F_struc(coneRows,1)'*Z(:);
lp.integer_variables = [];
lp.binary_variables = [];
lp.semicont_variables = [];
lp.x0 = [];
lp.options.saveduals = 1;
lp.options.savesolverinput = 0;
lp.options.savesolveroutput = 0;
lp.solver = p.solver.lower;
linearResult = feval(p.solver.lower.call,lp);
info.lpTime = linearResult.solvertime;
boundRows = nnz(isfinite(p.ub))+nnz(isfinite(p.lb));
linearDualRows = boundRows+(1:p.K.l);
if linearResult.problem~=0 || ~isfield(linearResult,'Dual') || ...
        numel(linearResult.Dual)<boundRows+p.K.l
    info.reason = 'original LP has no dual';
    return
end
At = -p.F_struc(:,2:end);
c = p.F_struc(:,1);
[bound,certificate] = bnbCertifySDPDualLPBound(At,-p.c,c,p.K,B, ...
    linearResult.Dual(linearDualRows),[],[],p.lb,p.ub);
% The certificate omits the constant objective term. Leave a rounding
% allowance when restoring it before comparing against the incumbent.
bound = bound+p.f-8*eps(1)*max(1,abs(p.f));
info.reason = certificate.reason;
if isfinite(bound)
    prune = bound >= upper+1e-6*max(1,abs(upper));
    if ~prune
        candidate = compressed.Primal;
        for g = 1:numel(groups)
            members = groups{g}(:);
            level = min(max(candidate(members(1)),0),numel(members));
            candidate(members) = min(max(level-(0:numel(members)-1)',0),1);
        end
        discrete = unique([p.binary_variables(:);p.integral_variables(:)]);
        if any(~isfinite(candidate)) || ...
                any(candidate<p.lb-1e-7) || any(candidate>p.ub+1e-7) || ...
                all(abs(candidate(discrete)- ...
                round(candidate(discrete)))<=p.options.bnb.inttol)
            % Integer candidates need an exact SDP node solve before BNB
            % can accept an incumbent or fathom the node.
            candidate = [];
        end
    end
end
end
