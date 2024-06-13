function varargout = cabs(varargin)

switch class(varargin{1})

    case 'double'
        varargout{1} = abs(varargin{1});        

    case 'sdpvar'
        varargout{1} = InstantiateElementWise(mfilename,varargin{:});

    case 'char'

        operator = CreateBasicOperator('convex','positive','callback');
        operator.derivative = @derivative;        
        operator.bounds = @bounds;                        
        operator.convexhull = @convexhull;  
        
        varargout{1} = [];
        varargout{2} = operator;
        varargout{3} = varargin{3};

    otherwise
        error([upper(mfilename) ' called with weird argument']);
end

function df = derivative(x)
df = zeros(length(x),1);
df(x>0) = 1;
df(x<0) = -1;

% Bounding functions for the branch&bound solver
function [L,U] = bounds(xL,xU)

if xL <= 0 && xU >= 0
    L = 0;
    U = max(abs(xL),abs(xU));
elseif xL >= 0
    L = xL;
    U = xU;
else
    L = abs(xU);
    U = abs(xL);
end