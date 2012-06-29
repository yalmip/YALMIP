function F = linearnegativeconstraint_iff_binary(f,X,M,m,eps);
% Big-M for f<=0 iff X==1. Assumes f and X vectors of same size

if nargin < 3 | isempty(M)
    [M,m,infbound] = derivebounds(f);
    if infbound
        warning('You have unbounded variables in IFF leading to a lousy big-M relaxation.');
    end
end

if nargin < 5
    eps = 1e-5;
end

if length(X) ~= length(f)
    error('Inconsistent sizes in linearnegativeconstraint_implies_binary: Report bug');
end

% X == 1 implies f<=0
F = [f <= M.*(1-X)];

% X == 0 implies f>=0
F = [F, f >= m.*X];

% f < -eps implies X==1
F = [F, f >= -eps + (m+eps).*X];

% f > eps implies X == 0
F = [F, f <= eps+(M-eps).*(1-X)];