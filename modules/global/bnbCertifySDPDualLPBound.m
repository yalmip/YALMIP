function [bound,details] = bnbCertifySDPDualLPBound(At,b,c,K,B,y,extraAt,extraC,explicitLB,explicitUB)
%BNBCERTIFYSDPDUALLPBOUND Conservative SDP dual Lagrangian lower bound.
% B is a floating-point factor. The exact-real matrix B'*B is PSD even if
% the rounded product is not. An error interval covers each scalar product.
% The LP multipliers y may be approximate: after clipping to nonnegative,
% their stationarity residual is minimized explicitly over variable boxes.
% Returns -Inf when a finite bound cannot be established.
if nargin < 7, extraAt = sparse(0,numel(b)); extraC = zeros(0,1); end
if nargin < 9
    explicitLB = -inf(numel(b),1);
    explicitUB = inf(numel(b),1);
end
At = [At(1:K.l,:);extraAt;At(K.l+1:end,:)];
c = [c(1:K.l);extraC;c(K.l+1:end)];
m = K.l+size(extraAt,1);
n = numel(b);
linear = At(1:m,:);
linearC = c(1:m);
cone = At(m+1:end,:);
coneC = c(m+1:end);
bound = -inf;
details = struct('reason','invalid input','maxResidual',NaN, ...
    'roundingBudget',NaN);
if numel(K.s)~=1 || size(B,2)~=K.s(1) || ...
        size(At,1)~=m+K.s(1)^2 || size(At,2)~=n || ...
        numel(c)~=size(At,1) || numel(explicitLB)~=n || ...
        numel(explicitUB)~=n || ~isreal(At) || ~isreal(b) || ...
        ~isreal(c) || ~isreal(B) || ~isreal(y) || ...
        any(isnan(explicitLB(:))) || any(isnan(explicitUB(:))) || ...
        numel(y)~=m || any(~isfinite(B(:))) || any(~isfinite(y(:))) || ...
        any(~isfinite(nonzeros(At))) || any(~isfinite(c)) || ...
        any(~isfinite(b))
    return
end
dimension = size(B,2);
unit = eps(1);
gamma = @(k) (k*unit)./(1-k*unit);
if max([size(B,1),dimension,numel(coneC),m,n])*unit >= 0.01
    return
end
% Keep every nonzero product well inside the normal double range. Outside
% this range, an underflowed dot product needs a separate absolute error
% analysis; an exact node solve is safer and inexpensive as a fallback.
finiteBox = [explicitLB(isfinite(explicitLB)); ...
    explicitUB(isfinite(explicitUB))];
inputScale = abs([nonzeros(At);c(:);b(:);B(:);y(:);finiteBox(:)]);
inputScale = inputScale(inputScale>0);
if any(inputScale<1e-50 | inputScale>1e50)
    details.reason = 'extreme input scale';
    return
end
% Construct a PSD dual as an exact Gram product of stored doubles.
Z = B'*B;
Zerror = 4*gamma(size(B,1))*(abs(B)'*abs(B));
if any(~isfinite(Z(:))) || any(~isfinite(Zerror(:)))
    details.reason = 'nonfinite Gram product';
    return
end
z = Z(:);
ze = Zerror(:);
y = max(0,y(:));
% A one-variable row gives an outward box bound. The box need not be tight.
lb = -inf(n,1);
ub = inf(n,1);
for row = 1:m
    [~,column,coefficient] = find(linear(row,:));
    if numel(column)~=1, continue; end
    ratio = linearC(row)/coefficient;
    radius = 8*unit*max(1,abs(ratio));
    if coefficient > 0
        ub(column) = min(ub(column),ratio+radius);
    else
        lb(column) = max(lb(column),ratio-radius);
    end
end
lb = max(lb,explicitLB(:)-8*unit*max(1,abs(explicitLB(:))));
ub = min(ub,explicitUB(:)+8*unit*max(1,abs(explicitUB(:))));
if any(isnan(lb)) || any(isnan(ub)) || any(lb>ub)
    details.reason = 'inconsistent box';
    return
end
% Affine Lagrangian f(x)-<Z,S(x)>-y'*(c_linear-A_linear*x).
q = -b(:) + cone'*z;
qAbsolute = abs(cone)'*abs(z);
qError = abs(cone)'*ze + ...
    8*gamma(full(sum(spones(cone),1))') .* qAbsolute + ...
    8*unit*(abs(b(:))+abs(q));
% Interior-point LP duals often have a tiny reduced cost in a one-sided
% unbounded direction. Decrease the multiplier of that variable's explicit
% one-variable bound row to move the whole interval to the safe side.
for pass = 1:2
    provisional = q+linear'*y;
    provisionalError = qError + ...
        8*gamma(full(sum(spones(linear),1))') .* ...
        (abs(linear)'*abs(y)) + ...
        8*unit*(abs(q)+abs(provisional));
    changed = false;
    for j = 1:n
        if ~isfinite(ub(j)) && provisional(j)-provisionalError(j)<0
            rows = find(linear(:,j)<0 & ...
                full(sum(spones(linear),2))==1 & y>0);
            if isempty(rows), continue; end
            row = rows(1);
            delta = (4*provisionalError(j)-provisional(j))/ ...
                abs(full(linear(row,j)));
            if delta>0 && delta<y(row)
                y(row) = y(row)-delta;
                changed = true;
            end
        elseif ~isfinite(lb(j)) && provisional(j)+provisionalError(j)>0
            rows = find(linear(:,j)>0 & ...
                full(sum(spones(linear),2))==1 & y>0);
            if isempty(rows), continue; end
            row = rows(1);
            delta = (4*provisionalError(j)+provisional(j))/ ...
                full(linear(row,j));
            if delta>0 && delta<y(row)
                y(row) = y(row)-delta;
                changed = true;
            end
        end
    end
    if ~changed, break; end
end
constant = -coneC'*z - linearC'*y;
constantAbsolute = abs(coneC)'*abs(z) + abs(linearC)'*abs(y);
constantError = abs(coneC)'*ze + ...
    8*gamma(nnz(coneC)+m)*constantAbsolute + ...
    8*unit*abs(constant);
reduced = q+linear'*y;
reducedAbsolute = abs(linear)'*abs(y);
reducedError = qError + ...
    8*gamma(full(sum(spones(linear),1))') .* reducedAbsolute + ...
    8*unit*(abs(q)+abs(reduced));
if any(~isfinite(reducedError)) || ~isfinite(constantError)
    details.reason = 'nonfinite interval';
    return
end
lowerReduced = reduced-reducedError;
upperReduced = reduced+reducedError;
% For an unbounded direction, the entire coefficient interval must have
% the favorable sign. Otherwise the certificate cannot prove a bound.
if any((~isfinite(ub) & lowerReduced<0) | ...
       (~isfinite(lb) & upperReduced>0) | ...
       (~isfinite(lb) & ~isfinite(ub)))
    details.reason = 'unbounded residual direction';
    details.maxResidual = max(abs(reduced));
    details.unboundedIndex = find(~isfinite(ub) | ~isfinite(lb));
    details.unboundedResidual = reduced(details.unboundedIndex);
    details.unboundedError = reducedError(details.unboundedIndex);
    return
end
contribution = zeros(n,1);
for j = 1:n
    choices = [];
    if isfinite(lb(j))
        choices = [lowerReduced(j)*lb(j),upperReduced(j)*lb(j)];
    end
    if isfinite(ub(j))
        choices = [choices,lowerReduced(j)*ub(j),upperReduced(j)*ub(j)]; %#ok<AGROW>
    end
    if isempty(choices)
        details.reason = 'missing finite endpoint';
        return
    end
    contribution(j) = min(choices);
end
roundingBudget = 8*gamma(n+1)* ...
    (abs(constant)+constantError+sum(abs(contribution)));
bound = constant-constantError+sum(contribution)-roundingBudget;
if ~isfinite(bound)
    bound = -inf;
    details.reason = 'nonfinite bound';
    return
end
details.reason = 'finite';
details.maxResidual = max(abs(reduced));
details.roundingBudget = roundingBudget+constantError;
end
