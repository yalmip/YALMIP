function F = binary_implies_linearequality(f,X,M,m,eps);

if nargin < 3 || isempty(M)
    [M,m,infbound] = derivebounds(f);
    if infbound
        warning('You have unbounded variables in an implication leading to a lousy big-M relaxation.');
    end
end

if all(m==0)
    lhs = m;
elseif all(m==1)
    lhs = 1-X;
else
    lhs = m.*(1-X);
end
if all(M==0)
    rhs = M;
elseif all(M==1)
    rhs = 1-X;
else
    rhs = M.*(1-X);
end
F = [f -f] <= [rhs -lhs];