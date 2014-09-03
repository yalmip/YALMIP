function varargout = or(varargin)
%OR (overloaded)

% Models OR using a nonlinear operator definition
if nargin == 1
    varargout{1} = varargin{1};
else
    varargout{1} = (yalmip('define','lmior',varargin{:}) == 1);
end
