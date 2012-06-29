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
% The model for OR is set(z>x) + set(z>y) + set(z<x+y) + set(binary(z))
%
% It is assumed that x and y are binary variables (either explicitely
% declared using BINVAR, or constrained using BINARY.)
%
%   See also SDPVAR/AND, BINVAR, BINARY

% Author Johan Löfberg 
% $Id: or.m,v 1.1 2006-08-10 18:00:21 joloef Exp $   

% Author Johan Löfberg 
% $Id: or.m,v 1.1 2006-08-10 18:00:21 joloef Exp $   

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
        
        varargout{1} = set(sum(xy) > z) + set(z > xy) +set(binary(z)) ;
        varargout{2} = struct('convexity','milp','monotonicity','milp','definiteness','milp');
        varargout{3} = xy;

    case 'sdpvar'
        x = varargin{1};
        y = varargin{2};
        varargout{1} = yalmip('addextendedvariable','or',varargin{:});

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