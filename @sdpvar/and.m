function varargout = and(varargin)
%AND (overloaded)
%   
%    z = and(x,y)
%    z = x & y
%
% It is assumed that x and y are binary variables (either explicitely
% declared using BINVAR, or constrained using BINARY.)
%
% See also SDPVAR/OR, SDPVAR/XOR, SDPVAR/NOT, BINVAR, BINARY

if isa(varargin{1},'double')
    varargin = {varargin{[2:nargin 1]}};
end
switch class(varargin{1})
    case 'char'
        
        z = varargin{2};        
        xy = [];        
        for i = 3:nargin            
            xy = [xy varargin{i}];
        end

        varargout{1} = (xy >= z) + (length(xy)-1+z >= sum(xy)) + (binary(z));
        varargout{2} = struct('convexity','none','monotonicity','none','definiteness','none','model','integer');
        varargout{3} = xy;

    case 'sdpvar'
        varargout{1} = yalmip('define','and',varargin{:});

    otherwise
end