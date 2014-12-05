function varargout = not(varargin)
%NOT (overloaded)
%   
%    z = ~x
%    z = not(x)
%
% It is assumed that x is a binary variables (either explicitely
% declared using BINVAR, or constrained using BINARY.)
%
% See also SDPVAR/AND, SDPVAR/OR, SDPVAR/XOR, BINVAR, BINARY

switch class(varargin{1})
    case 'char'
        
        z = varargin{2};
        x = varargin{3};
        
        varargout{1} = z == 1-x;
        varargout{2} = struct('convexity','none','monotonicity','none','definiteness','none','model','integer');
        varargout{3} = recover(getvariables(x));
        
    case 'sdpvar'
        if nargin > 1
            error('Wrong numbers of arguments.');
        end
        varargout{1} = yalmip('define','not',varargin{1});
        
    otherwise
end