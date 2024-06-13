function [F,z_normalizing,X] = NormalizeCallback(varargin)

X = [];
doAssignInitials = varargin{end};
z_normalizing = varargin{end-1};
for i = 3:nargin-2
    if isa(varargin{i},'sdpvar')
        X = varargin{i};
        break
    end
end
n = length(X);
if isequal(getbase(X),[spalloc(n,1,0) speye(n)])
    z_normalizing = [];
    F = lmi([]);
else
    if doAssignInitials
        dX = value(X);
        if ~all(isnan(dX))
            assign(z_normalizing,dX);
        end
    end
    try        
        F = X == z_normalizing;
    catch
        disp('Report bug in NORMALIZECALLBACK');
        error('Report bug in NORMALIZECALLBACK')
    end
end
