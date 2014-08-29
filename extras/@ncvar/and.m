function varargout = and(varargin)
%AND (overloaded)
%   
%    z = and(x,y)
%    z = x & y
%
% The AND operator is implemented using the concept of nonlinear operators
% in YALMIP. X|Y defines a new so called derived variable that can be
% treated as any other variable in YALMIP. When SOLVESDP is issued,
% constraints are added to the problem to model the AND operator. The new
% constraints add constraints to ensure that z, x and y satisfy the
% truth-table for AND.
%
% The model for AND is (z<=x)+(z<=y)+(1+z>=x+y)+(binary(z))
%
% It is assumed that x and y are binary variables (either explicitely
% declared using BINVAR, or constrained using BINARY.)
%
%   See also SDPVAR/AND, BINVAR, BINARY

switch class(varargin{1})
    case 'char'
        z = varargin{2};
        x = varargin{3};
        y = varargin{4};
    
        % *******************************************************
        % For *some* efficiency,we merge expressions like A&B&C&D
        xvars = getvariables(x);
        yvars = getvariables(y);
        allextvars = yalmip('extvariables');
        
        if (length(xvars)==1) & ismembc(xvars,allextvars)
            x = expandand(x,allextvars);
        end
        
        if (length(yvars)==1) & ismembc(yvars,allextvars)
            y = expandand(y,allextvars);
        end
        % *******************************************************
                
        varargout{1} = (x >= z) + (y >= z) + (length(x)+length(y)-1+z >= sum(x)+sum(y)) + (binary(z));
        varargout{2} = struct('convexity','milp','monotoncity','milp','definiteness','milp');
        varargout{3} = [];

    case 'sdpvar'
        x = varargin{1};
        y = varargin{2};

        varargout{1} = yalmip('addextendedvariable','and',varargin{:});

    otherwise
end

function x = expandand(x,allextvars)

xmodel = yalmip('extstruct',getvariables(x));

if isequal(xmodel.fcn,'and')
    x1 = xmodel.arg{1};
    x2 = xmodel.arg{2};
    if  ismembc(getvariables(xmodel.arg{1}),allextvars)
        x1 = expandand(xmodel.arg{1},allextvars);     
    end
    if  ismembc(getvariables(xmodel.arg{2}),allextvars)
        x2 = expandand(xmodel.arg{2},allextvars);     
    end
    x = [x1 x2];
end