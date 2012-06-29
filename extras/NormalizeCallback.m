function F = NormalizeCallback(varargin)

z_normalizing = varargin{end};
X = varargin{3};
n = length(X);
if isequal(getbase(X),[zeros(n,1) eye(n)])
    F = set([]);
else
    dX = double(X);
    if ~all(isnan(dX))
        assign(z_normalizing,double(X));
    end
    try
        %[M,m] = derivebounds(X);
        F = [X == z_normalizing];
    catch
        disp('Report bug in NORMALIZECALLBACK');
        error('Report bug in NORMALIZECALLBACK')
    end
end