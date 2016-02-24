function F = NormalizeCallback(varargin)

z_normalizing = varargin{end};
for i = 3:nargin-1
    if isa(varargin{i},'sdpvar')
        X = varargin{i};
        break
    end
end
%X = varargin{3};
n = length(X);
if isequal(getbase(X),[spalloc(n,1,0) speye(n)])
    F = lmi([]);
else
    dX = value(X);
    if ~all(isnan(dX))
        assign(z_normalizing,dX);
    end
    try
        %[M,m] = derivebounds(X);
        F = X == z_normalizing;
    catch
        disp('Report bug in NORMALIZECALLBACK');
        error('Report bug in NORMALIZECALLBACK')
    end
end
