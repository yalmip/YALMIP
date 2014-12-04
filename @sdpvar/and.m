function varargout = and(varargin)
%ANd (overloaded)
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
% The model for ARE is (z<=x)+(z<=y)+(1+z>=x+y)+(binary(z))
%
% It is assumed that x and y are binary variables (either explicitely
% declared using BINVAR, or constrained using BINARY.)
%
%   See also SDPVAR/AND, BINVAR, BINARY

if isa(varargin{1},'double')
    varargin = {varargin{[2:nargin 1]}};
end
switch class(varargin{1})
    case 'char'
        
        z = varargin{2};        
        xy = [];
        allextvars = yalmip('extvariables');
        for i = 3:nargin
            x = varargin{i};
            xvars = getvariables(x);
            if (length(xvars)==1) & ismembc(xvars,allextvars)
                x = expandand(x,allextvars);
            end
            xy = [xy x];
        end

        varargout{1} = (xy >= z) + (length(xy)-1+z >= sum(xy)) + (binary(z));
        varargout{2} = struct('convexity','none','monotonicity','none','definiteness','none','model','integer');
        varargout{3} = xy;

    case 'sdpvar'
        varargout{1} = yalmip('define','and',varargin{:});

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