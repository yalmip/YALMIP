function varargout = or(varargin)
%OR (overloaded)
%   
%    z = or(x,y)
%    z = x | y
%
% The OR operator is implemented using the concept of nonlinear operators
% in YALMIP. X|Y defines a new so called derived variable that can be
% treated as any other variable in YALMIP. When SOLVESDP is issued,
% constraints are added to the problem to model the OR operator. The new
% constraints add constraints to ensure that z,x and y satisfy the
% truth-table for OR. 
%
% The model for OR is [z>=x, z>=y, z<=x+y, binary(z)]
%
% It is assumed that x and y are binary variables (either explicitely
% declared using BINVAR, or constrained using BINARY.)
%
%   See also SDPVAR/AND, BINVAR, BINARY

% Author Johan Löfberg 
% $Id: or.m,v 1.15 2007-08-02 19:17:36 joloef Exp $   

% Author Johan Löfberg 
% $Id: or.m,v 1.15 2007-08-02 19:17:36 joloef Exp $   

% Models OR using a nonlinear operator definition
switch class(varargin{1})
    case 'char'
        z = varargin{2};
        x = varargin{3};
        y = varargin{4};        
          
        % *******************************************************
        % For *some* efficiency,we merge expressions like A|B|C|D
        xvars = getvariables(x);
        yvars = getvariables(y);
        allextvars = yalmip('extvariables');
        
        if (length(xvars)==1) & ismembc(xvars,allextvars)
            x = expandor(x,allextvars);
        end
        
        if (length(yvars)==1) & ismembc(yvars,allextvars)
            y = expandor(y,allextvars);
        end
        % *******************************************************
                 
        xy=[x y];
        
        varargout{1} = set(sum(xy) >= z) + set(z >= xy) +set(binary(z));
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
            x = varargin{1};
            y = varargin{2};
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