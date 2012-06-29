function varargout = sqrtm(varargin)
%SQRTM (overloaded)

% Author Johan Löfberg
% $Id: sqrtm.m,v 1.2 2006-11-15 10:03:57 joloef Exp $

switch nargout
    case 0
        sqrtm_internal(varargin{:});
    case 1
        varargout{1} = sqrtm_internal(varargin{:});
    case 2
        [varargout{1},varargout{2}] = sqrtm_internal(varargin{:});
    case 3
        [varargout{1},varargout{2},varargout{3}] = sqrtm_internal(varargin{:});
    otherwise
end
