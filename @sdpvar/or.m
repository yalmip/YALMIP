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

    case 'sdpvar'
        if nargin == 1
            if length(varargin{1})==1
                varargout{1} = varargin{1}
            else
                x = varargin{1};
                % bug in matlab...
                %temp = or(x(1),x(2));
                temp = or(extsubsref(x,1),extsubsref(x,2));
                for i = 3:length(x)
                    temp = or(temp,extsubsref(x,i));
                end
                 varargout{1} = temp;
            end           
        else          
            varargout{1} = yalmip('define','or',varargin{:});
        end
    otherwise
end