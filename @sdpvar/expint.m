function varargout = expint(varargin)
%EXPINT (overloaded)

% Author Johan Löfberg
switch class(varargin{1})

    case 'sdpvar'
        varargout{1} = InstantiateElementWise(mfilename,varargin{:});

    case 'char'
        
        varargout{1} = [];
        varargout{2} = createOperator;
        varargout{3} = varargin{3};

    otherwise
        error('SDPVAR/EXPINT called with CHAR argument?');
end

function operator = createOperator

operator = struct('convexity','convex','monotonicity','decreasing','definiteness','positive','model','callback');
operator.bounds     = @bounds;
operator.derivative = @(x)(-exp(-x)./x);
operator.range = [0 inf];
operator.domain = [1e-8 inf];

% Bounding functions for the branch&bound solver
function [L,U] = bounds(xL,xU)
L = expint(xU);
U = expint(xL);
