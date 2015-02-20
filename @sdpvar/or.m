function varargout = or(varargin)
%OR (overloaded)
%   
%    z = or(x,y)
%    z = x | y
%
% It is assumed that x and y are binary variables (either explicitely
% declared using BINVAR, or constrained using BINARY.)
%
% See also SDPVAR/AND, SDPVAR/XOR, SDPVAR/NOT, BINVAR, BINARY

switch class(varargin{1})
    case 'char'
        
        z = varargin{2};        
        xy = [];        
        for i = 3:nargin            
            xy = [xy varargin{i}];
        end
                
        varargout{1} = (sum(xy) >= z) + (z >= xy) +(binary(z));
        varargout{2} = struct('convexity','none','monotonicity','none','definiteness','none','model','integer');
        varargout{3} = xy;
            
    case {'sdpvar','double','logical'}
        
        varargout{1} = vectorizedlogic(@or,varargin{:});

    otherwise
end