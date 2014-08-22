function F = binary_implies_linearnegativeconstraint(f,X,M,m,eps);
% Big-M for f<=0 iff X==1. Assumes f and X vectors of same size

if nargin < 3 || isempty(M)
    [M,m,infbound] = derivebounds(f);
    if infbound
        warning('You have unbounded variables in an implication leading to a lousy big-M relaxation.');
    end
end

if nargin < 5
    eps = 1e-5;
end

F = [f <= M.*(1-X)];