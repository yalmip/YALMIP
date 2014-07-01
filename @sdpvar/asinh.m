function varargout = asinh(varargin)
%ASINH (overloaded)

switch class(varargin{1})

    case 'sdpvar'
        varargout{1} = InstantiateElementWise(mfilename,varargin{:});

    case 'char'

        operator = struct('convexity','none','monotonicity','increasing','definiteness','none','model','callback');
        operator.convexhull = [];
        operator.bounds = @bounds;
        operator.derivative = @(x)((1 + x.^2).^-0.5);
            
        varargout{1} = [];
        varargout{2} = operator;
        varargout{3} = varargin{3};

    otherwise
        error('SDPVAR/ASINH called with CHAR argument?');
end

function [L,U] = bounds(xL,xU)
L = asinh(xL);
U = asinh(xU);