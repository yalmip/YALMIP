function varargout = acosh(varargin)
%ACOSH (overloaded)

switch class(varargin{1})

    case 'sdpvar'
        varargout{1} = InstantiateElementWise(mfilename,varargin{:});

    case 'char'

        operator = struct('convexity','none','monotonicity','none','definiteness','positive','model','callback');
        operator.convexhull = [];
        operator.bounds = @bounds;

        varargout{1} = [];
        varargout{2} = operator;
        varargout{3} = varargin{3};

    otherwise
        error('SDPVAR/ACOSH called with CHAR argument?');
end

function [L,U] = bounds(xL,xU)
if xU<=-1
    L = acosh(xU);
    U = acosh(xL);
elseif xL>=1
    L = acosh(xL);
    U = acosh(xU);
else
    L = -inf;
    U = inf;
end

