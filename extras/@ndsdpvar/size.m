function varargout = size(X,d)
% SIZE (overloaded)

if nargin == 1 && nargout > 1 && nargout == length(X.dim)
    for i = 1:length(X.dim);
        varargout{i} = X.dim(i);
    end
    return
end

if nargin == 1
    varargout{1} = X.dim;
else
    if (d < 1) || d~=round(d)
        error('Dimension argument must be a positive integer scalar in the range 1 to 2^31.');
    end
    if d <= length(X.dim)
        varargout{1} = X.dim(d);
    else
        varargout{1} = 1;
    end
end