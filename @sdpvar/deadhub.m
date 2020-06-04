function varargout = deadhub(varargin)

switch class(varargin{1})

    case 'sdpvar'
        varargout{1} = InstantiateElementWise(mfilename,varargin{:});

    case 'char'

        operator = CreateBasicOperator('convex','callback');
        operator.bounds     = @bounds;
        operator.convexhull = @convexhull;      

        varargout{1} = [];
        varargout{2} = operator;
        varargout{3} = varargin{3};

    otherwise
        error(['SDPVAR/' upper(mfilename) ' called with weird argument']);
end

function [L,U] = bounds(xL,xU,lambda)

fL = deadhub(xL,lambda);
fU = deadhub(xU,lambda);
U = max(fL,fU);
L = min(fL,fU);
if xL<0 & xU>0
    L = 0;
end

function [Ax, Ay, b, K] = convexhull(xL,xU,lambda)

K.l = 0;
K.f = 0;
fL = deadhub(xL,lambda);
fU = deadhub(xU,lambda);
if xL>=-lambda & xU<=lambda
     Ax = 0;Ay = 1;b = 0;K.f = 1;
elseif xU < -3*lambda
    Ax = 1;Ay = 1;b = 2*lambda^2;K.f = 1;
elseif xL > 3*lambda
    Ax = -1;Ay = 1;b = 2*lambda^2;K.f = 1;
else
    dfL = derivative(xL,lambda);
    dfU = derivative(xU,lambda);
    [Ax,Ay,b,K] = convexhullConvex(xL,xU,fL,fU,dfL,dfU);
end

function df=derivative(x,lambda)
ax = abs(x);
if ax<lambda
    df=0;
elseif ax>3*lambda
    df = lambda;
elseif ax<=3*lambda
    df = 0.25*(2*ax-6*lambda)+lambda;
end
if x<0
    df=-df;
end