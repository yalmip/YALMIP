function s = size(X,d)
% SIZE (overloaded)

if nargin == 1
    s = X.dim;
else
    if (d < 1) || d~=round(d)
        error('Dimension argument must be a positive integer scalar in the range 1 to 2^31.');
    end
    if d <= length(X.dim)
        s = X.dim(d);
    else
        s = 1;
    end
end