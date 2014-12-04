function varargout = not(varargin)
%NOT (overloaded)
%   
%    z = ~x
%    z = not(x)
%
% It is assumed that x and y are binary variables (either explicitely
% declared using BINVAR, or constrained using BINARY.)
%
% See also SDPVAR/ANS, SDPVAR/OR, SDPVAR/XOR, BINVAR, BINARY

switch class(varargin{1})
    case 'char'
        
        z = varargin{2};        
        xy = [];
        allextvars = yalmip('extvariables');
     
        x = varargin{3};
        xvars = getvariables(x);
        if (length(xvars)==1) & ismembc(xvars,allextvars)
            x = expandnot(x,allextvars);
        end
                
        varargout{1} = z == 1-x;
        varargout{2} = struct('convexity','none','monotonicity','none','definiteness','none','model','integer');
        varargout{3} = x;

    case 'sdpvar'
        if nargin > 1
            error('Wrong numbers of arguments.');
        end
        varargout{1} = yalmip('define','not',varargin{1});

    otherwise
end

function x = expandnot(x,allextvars)

xmodel = yalmip('extstruct',getvariables(x));

if isequal(xmodel.fcn,'not')
    x = [];
    for i = 1:length(xmodel.arg)
        xi = xmodel.arg{i};
        if  ismembc(getvariables(xi),allextvars)
            xi = expandnot(xi,allextvars);
        end
        x = [x xi];
    end
end