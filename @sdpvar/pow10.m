function varargout = pow10(varargin)
%POW10 (overloaded)

switch class(varargin{1})

    case 'double'
        error('Overloaded SDPVAR/POW10 CALLED WITH DOUBLE. Report error')

    case 'sdpvar'
        varargout{1} = InstantiateElementWise(mfilename,varargin{:});

    case 'char'

        operator = struct('convexity','convex','monotonicity','increasing','definiteness','positive','model','callback');
        operator.convexhull = @convexhull;
        operator.bounds     = @bounds;
        operator.derivative = @derivative;

        varargout{1} = [];
        varargout{2} = operator;
        varargout{3} = varargin{3};

    otherwise
        error('SDPVAR/POW10 called with CHAR argument?');
end

function df = derivative(x)
df = log(10)*pow10(x);

function [L,U] = bounds(xL,xU)
L = pow10(xL);
U = pow10(xU);

function [Ax, Ay, b] = convexhull(xL,xU)
fL = pow10(xL);
fU = pow10(xU);
dfL = log(10)*fL;
dfU = log(10)*fU;
[Ax,Ay,b] = convexhullConvex(xL,xU,fL,fU,dfL,dfU);