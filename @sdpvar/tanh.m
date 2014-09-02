function varargout = tanh(varargin)
%TANH (overloaded)

switch class(varargin{1})

    case 'double'
        error('Overloaded SDPVAR/TANH CALLED WITH DOUBLE. Report error')

    case 'sdpvar'
        varargout{1} = InstantiateElementWise(mfilename,varargin{:});

    case 'char'

        operator = struct('convexity','none','monotonicity','increasing','definiteness','none','model','callback');
        operator.convexhull = [];
        operator.bounds = @bounds;
        operator.derivative = @(x)(1-tanh(x).^2);
        operator.range = [-1 1];

        varargout{1} = [];
        varargout{2} = operator;
        varargout{3} = varargin{3};

    otherwise
        error('SDPVAR/TANH called with CHAR argument?');
end

function [L,U] = bounds(xL,xU)
L = cosh(xL);
U = cosh(xU);

