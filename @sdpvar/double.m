function [sys,values] = double(varargin)
% OBSOLETE. Use VALUE instead

% New syntax
switch nargout
    case 0
        value(varargin{:})
    case 1
        sys = value(varargin{:});
    case 2
        [sys,values] = value(varargin{:});
    otherwise
        error('Too many output arguments.');
end

