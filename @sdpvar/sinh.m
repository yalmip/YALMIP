function varargout = sinh(varargin)
%SINH (overloaded)

switch class(varargin{1})

    case 'double'
        error('Overloaded SDPVAR/ACOT CALLED WITH DOUBLE. Report error')

    case 'sdpvar'
        varargout{1} = InstantiateElementWise(mfilename,varargin{:});

    case 'char'

        operator = struct('convexity','none','monotonicity','none','definiteness','none','model','callback');
        operator.convexhull = [];
        operator.bounds = @bounds;
        operator.derivative = @(x)(cosh(x));

        varargout{1} = [];
        varargout{2} = operator;
        varargout{3} = varargin{3};

    otherwise
        error('SDPVAR/SINH called with CHAR argument?');
end

function [L,U] = bounds(xL,xU)
L = sinh(xL);
U = sinh(xU);
