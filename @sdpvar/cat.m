function y = cat(varargin)
%CAT (overloaded)

switch varargin{1}
    case 1
        y = vertcat(varargin{2:end});
    case 2
        y = horzcat(varargin{2:end});
    otherwise
        for i = 2:nargin
            if isa(varargin{i},'sdpvar');
                varargin{i} = ndsdpvar(varargin{i});
            end
        end
        y = cat(varargin{:});
end