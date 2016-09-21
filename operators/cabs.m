function varargout = cabs(varargin)
%cabs (overloaded)

switch class(varargin{1})

    case 'double'
        varargout{1} = abs(varargin{1});        

    case 'sdpvar'
        varargout{1} = InstantiateElementWise(mfilename,varargin{:});

    case 'char'

        operator = struct('convexity','convex','monotonicity','none','definiteness','positive','model','callback');
        operator.derivative = @derivative;        
        operator.bounds = @bounds;                        
        operator.convexhull = @convexhull;                        
        varargout{1} = [];
        varargout{2} = operator;
        varargout{3} = varargin{3};

    otherwise
        error('SDPVAR/CABS called with CHAR argument?');
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

function [Ax, Ay, b, K] = convexhull(xL,xU)
fL = abs(xL);
fU = abs(xU);
if fL == fU || xL>=0
    Ax = -1;
    Ay = 1;
    b = 0;
    K.f = 1;
    K.l = 0;
elseif xU<=0
    Ax = 1;
    Ay = 1;
    b = 0;
    K.f = 1;
    K.l = 0;
elseif xL < 0 && xU > 0
    dfL = -1;
    dfU = 1;   
    [Ax,Ay,b] = convexhullConvex(xL,xU,fL,fU,dfL,dfU);
    K = [];
end
