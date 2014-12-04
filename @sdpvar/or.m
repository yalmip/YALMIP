function varargout = or(varargin)
%OR (overloaded)
%   
%    z = or(x,y)
%    z = x | y
%
% The model for OR is [z>=x, z>=y, z<=x+y, binary(z)]
%
% It is assumed that x and y are binary variables (either explicitely
% declared using BINVAR, or constrained using BINARY.)
%
% See also SDPVAR/AND, SDPVAR/XOR, SDPVAR/NOT, BINVAR, BINARY

switch class(varargin{1})
    case 'char'
        z = varargin{2};
        
        xy = [];
        allextvars = yalmip('extvariables');
        for i = 3:nargin
            x = varargin{i};
            xvars = getvariables(x);
            if (length(xvars)==1) & ismembc(xvars,allextvars)
                x = expandor(x,allextvars);
            end
            xy = [xy x];
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

function x = expandor(x,allextvars)

xmodel = yalmip('extstruct',getvariables(x));

if isequal(xmodel.fcn,'or')
    x1 = xmodel.arg{1};
    x2 = xmodel.arg{2};
    if  ismembc(getvariables(xmodel.arg{1}),allextvars)
        x1 = expandor(xmodel.arg{1},allextvars);     
    end
    if  ismembc(getvariables(xmodel.arg{2}),allextvars)
        x2 = expandor(xmodel.arg{2},allextvars);     
    end
    x = [x1 x2];
end