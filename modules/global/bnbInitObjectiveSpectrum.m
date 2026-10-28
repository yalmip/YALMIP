function [spectrum,info] = bnbInitObjectiveSpectrum(p,upper)
% Enumerate a superset of attainable discrete linear objective values.
% Nearby coefficients share a representative only with an explicit error
% interval covering the original coefficients. No model coefficient changes.
% Ignoring constraints enlarges the set of possible costs, which is safe
% for a lower bound. Cache only values below the first finite incumbent;
% bnbRoundObjectiveSpectrum caps every lift at that incumbent. All budgets
% are performance guards: exceeding one disables this optional tightening.
spectrum = [];
maxStates = 4096;
info = struct('reason','unsupported','entries',0,'classes',0,'seconds',0);
timer = tic;
if ~isfinite(upper) || nnz(p.Q) || ~isempty(p.evalMap) || ...
        any(p.variabletype) || ~isequal(p.K.m,0) || p.IntegerCost
    return
end
used = find(p.c);
if isempty(used) || ~all(ismember(used,p.integral_variables))
    return
end
c = full(p.c(used)); c = c(:);
lo = ceil(p.lb(used)); hi = floor(p.ub(used));
lo = lo(:); hi = hi(:);
if any(~isfinite([lo;hi;c])) || any(lo>hi) || ...
        any(abs([lo;hi])>1e6) || numel(c)>2000
    return
end
basePoint = lo;
basePoint(c<0) = hi(c<0);
base = p.f+c'*basePoint;
width = hi-lo;
weights = abs(c);
active = width>0;
weights = weights(active); width = width(active);
scale = 1+abs(p.f)+sum(abs(c).*max(abs(lo),abs(hi)));
rounding = 128*eps*max(1,numel(c))*scale;
if ~isfinite(scale) || rounding>1e-5*max(1,abs(upper))
    info.reason = 'rounding budget';
    return
end
[weights,order] = sort(weights); width = width(order);
representatives = []; counts = []; coefficientError = 0;
for i = 1:numel(weights)
    if isempty(representatives) || ...
            weights(i)-representatives(end)>1e-10*max(1,weights(i))
        representatives(end+1) = weights(i); %#ok<AGROW>
        counts(end+1) = width(i); %#ok<AGROW>
    else
        counts(end) = counts(end)+width(i);
        coefficientError = coefficientError+width(i)*abs(weights(i)-representatives(end));
    end
end
info.classes = numel(counts);
if numel(counts)>8
    info.reason = 'too many coefficient classes';
    return
end
radius = coefficientError+rounding;
cap = upper-base+radius;
values = 0;
for i = 1:numel(counts)
    maxCount = min(counts(i),max(0,floor((cap+rounding)/representatives(i))));
    if numel(values)*(maxCount+1)>min(1000000,16*maxStates)
        info.reason = 'enumeration budget';
        info.seconds = toc(timer);
        return
    end
    values = values(:)+(0:maxCount)*representatives(i);
    values = unique(values(values<=cap+rounding));
    if numel(values)>maxStates
        info.reason = 'spectrum budget';
        info.seconds = toc(timer);
        return
    end
end
spectrum = struct('values',base+values(:),'radius',radius,'cap',upper);
info.reason = 'active'; info.entries = numel(values); info.seconds = toc(timer);
end
