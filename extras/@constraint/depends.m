function LinearVariables = depends(varargin)

if nargin > 1
    LinearVariables = [];
    for i = 1:nargin
        LinearVariables_i = depends(varargin{i});
        LinearVariables = [LinearVariables;LinearVariables_i(:)];
    end
    LinearVariables = unique(LinearVariables);
    return
else
    LinearVariables = depends(lmi(varargin{1}));
end