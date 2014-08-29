function varargout = pow2(varargin)
%POW2 (overloaded)

switch class(varargin{1})

    case 'double'
        error('Overloaded SDPVAR/POW2 CALLED WITH DOUBLE. Report error')

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
        error('SDPVAR/POW2 called with CHAR argument?');
end

function df = derivative(x)
df = log(2)*pow2(x);

function [L,U] = bounds(xL,xU)
L = pow2(xL);
U = pow2(xU);

function [Ax, Ay, b] = convexhull(xL,xU)
fL = pow2(xL);
fU = pow2(xU);
dfL = log(2)*fL;
dfU = log(2)*fU;
[Ax,Ay,b] = convexhullConvex(xL,xU,fL,fU,dfL,dfU);