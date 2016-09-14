function sys = cat(varargin)

if nargin == 1
    sys = varargin{1};
    return
end

y = horzcat(varargin{:});