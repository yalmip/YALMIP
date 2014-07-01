function varargout = atan(varargin)
%ATAN (overloaded)

switch class(varargin{1})

    case 'sdpvar'
        varargout{1} = InstantiateElementWise(mfilename,varargin{:});

    case 'char'

        operator = struct('convexity','none','monotonicity','increasing','definiteness','none','model','callback');
        operator.convexhull = @convexhull;
        operator.bounds = @bounds;
        operator.derivative = @(x)((1+x.^2).^-1);
        operator.inverse = @(x)(tan(x));
        operator.range = [-pi/2 pi/2];
        
        varargout{1} = [];
        varargout{2} = operator;
        varargout{3} = varargin{3};

    otherwise
        error('SDPVAR/ATAN called with CHAR argument?');
end

function [L,U] = bounds(xL,xU)
L = atan(xL);
U = atan(xU);

function [Ax, Ay, b] = convexhull(xL,xU)
fL = atan(xL);
fU = atan(xU);
dfL = 1/(1+xL^2);
dfU = 1/(1+xU^2);
if xL >= 0
    % Concave region
    [Ax,Ay,b] = convexhullConcave(xL,xU,fL,fU,dfL,dfU);
elseif xU <= 0
    % Convex region
    [Ax,Ay,b] = convexhullConvex(xL,xU,fL,fU,dfL,dfU);
else
    % Changes convexity. We're lazy and let YALMIP sample instead
    Ax = [];
    Ay = [];
    b = [];
end