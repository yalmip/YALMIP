function ChanceConstraint = horzcat(varargin)

if nargin == 2
    ChanceConstraint = varargin{1};
    for i = 1:length(varargin{2}.P)
        ChanceConstraint.P{end+1} = varargin{2}.P{i};
        ChanceConstraint.level{end+1} = varargin{2}.level{i};
    end
else
    ChanceConstraint = horzcat(varargin{1},horzcat(varargin{2:end}));
end