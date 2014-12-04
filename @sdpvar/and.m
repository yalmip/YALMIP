function varargout = and(varargin)
%AND (overloaded)
%   
%    z = and(x,y)
%    z = x & y
%
% The model for AND is (z<=x)+(z<=y)+(1+z>=x+y)+(binary(z))
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
    x = [];
    for i = 1:length(xmodel.arg)
        xi = xmodel.arg{i};
        if  ismembc(getvariables(xi),allextvars)
            xi = expandand(xi,allextvars);
        end
        x = [x xi];
    end
end