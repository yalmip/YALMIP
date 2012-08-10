function varargout = or(varargin)
%OR (overloaded)

% Author Johan Löfberg 
% $Id: or.m,v 1.2 2007-08-02 19:33:16 joloef Exp $   

% Models OR using a nonlinear operator definition
if nargin == 1
    varargout{1} = varargin{1};
else
    varargout{1} = set(yalmip('define','lmior',varargin{:}) == 1);
end
