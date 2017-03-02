function varargout = logistic(varargin)
% LOGISTIC  Returns logistic function 1./(1+exp(-x))
%
% y = LOGISTIC(x)
%
% For a real vector x, LOGISTIC returns (1+exp(-x)).^-1

switch class(varargin{1})

    case 'double'
        x = varargin{1};
        varargout{1} = 1./(1+exp(-x));

    case 'sdpvar'
        varargout{1} = InstantiateElementWise(mfilename,varargin{:});

    case 'char'

        [M,m] = derivebounds(varargin{3});
        if m >= 0
            operator = struct('convexity','concave','monotonicity','increasing','definiteness','positive','model','callback');
        elseif M <= 0
            operator = struct('convexity','convex','monotonicity','increasing','definiteness','positive','model','callback');
        else
            operator = struct('convexity','none','monotonicity','increasing','definiteness','positive','model','callback');
        end
        operator.convexhull = @convexhull;
        operator.bounds     = @bounds;
        operator.derivative = @(x)logistic(x).*(1-logistic(x));
        operator.inverse    = @(x)(log(x)-log(1-x));
        operator.range = [0 1];

        varargout{1} = [];
        varargout{2} = operator;
        varargout{3} = varargin{3};

    otherwise
        error('SDPVAR/LOGISTIC called with CHAR argument?');
end

% Bounding functions for the branch&bound solver
function [L,U] = bounds(xL,xU)
L = 1./(1+exp(-xL));
U = 1./(1+exp(-xU));

function [Ax, Ay, b, K] = convexhull(xL,xU)
if xU <=0
    fL = logistic(xL);
    fU = logistic(xU);
    dfL = fL*(1-fL);
    dfU = fU*(1-fU);
    [Ax,Ay,b] = convexhullConvex(xL,xU,fL,fU,dfL,dfU);
elseif xL>=0
    fL = logistic(xL);
    fU = logistic(xU);
    dfL = fL*(1-fL);
    dfU = fU*(1-fU);
    [Ax,Ay,b] = convexhullConcave(xL,xU,fL,fU,dfL,dfU);
else
    Ax = [];
    Ay = [];
    b = [];
end
K = [];