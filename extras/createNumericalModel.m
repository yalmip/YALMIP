function p = createNumericalModel(varargin);

p = emptyNumericalModel;

if nargin >= 1
    p.F_struc = varargin{1};
end
if nargin >= 2
    p.K = varargin{2};
end
