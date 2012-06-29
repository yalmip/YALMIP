function varargout = exp(varargin)
%EXP (overloaded)

% Author Johan Löfberg
% $Id: exp.m,v 1.30 2009-03-11 09:45:32 joloef Exp $
switch class(varargin{1})

    case 'sdpvar'
        varargout{1} = InstantiateElementWise(mfilename,varargin{:});

    case 'char'
        
        varargout{1} = [];
        varargout{2} = createOperator;
        varargout{3} = varargin{3};

    otherwise
        error('SDPVAR/EXP called with CHAR argument?');
end

function operator = createOperator

operator = struct('convexity','convex','monotonicity','increasing','definiteness','positive','model','callback');
operator.convexhull = @convexhull;
operator.bounds     = @bounds;
operator.derivative = @(x)exp(x);
operator.range = [0 inf];

% Bounding functions for the branch&bound solver
function [L,U] = bounds(xL,xU)
L = exp(xL);
U = exp(xU);

function [Ax, Ay, b, K] = convexhull(xL,xU)
fL = exp(xL);
fU = exp(xU);
if fL == fU
    Ax = [];
    Ay = [];
    b = [];
else
    dfL = exp(xL);
    dfU = exp(xU);
    % A cut with tangent parallell to upper bound is very efficient
    xM = log((fU-fL)/(xU-xL));
    fM = exp(xM);
    dfM = exp(xM);
    [Ax,Ay,b] = convexhullConvex(xL,xM,xU,fL,fM,fU,dfL,dfM,dfU);
end
K = [];
