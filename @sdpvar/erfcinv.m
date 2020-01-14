function varargout = erfcinv(varargin)
%ERFCINV (overloaded)

switch class(varargin{1})

    case 'double'
        error('Overloaded SDPVAR/ERFCINV CALLED WITH DOUBLE. Report error')

    case 'sdpvar'
        varargout{1} = InstantiateElementWise(mfilename,varargin{:});

    case 'char'

        X = varargin{3};
        F = (-1+1e-9 <= X <= 1-1e-9);
        operator = struct('convexity','none','monotonicity','decreasing','definiteness','none','model','callback');
        operator.bounds = @bounds;
        operator.inverse = @(x)(erfc(x));
        operator.derivative = @(x)-1./(exp(-((erfcinv(x))).^2)*2/sqrt(pi))
        
        varargout{1} = [];
        varargout{2} = operator;
        varargout{3} = varargin{3};

    otherwise
        error('SDPVAR/ERFCINV called with CHAR argument?');
end

function [L,U] = bounds(xL,xU)
if xL<=-1
    U = inf
else
    U = erfcinv(xL);
end
if xU>=1
    L = -inf;
else
    L = erfcinv(xU);
end