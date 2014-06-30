function varargout = asec(varargin)
%ASEC (overloaded)

switch class(varargin{1})

    case 'sdpvar'
        varargout{1} = InstantiateElementWise(mfilename,varargin{:});

    case 'char'

        operator = struct('convexity','none','monotonicity','none','definiteness','none','model','callback');
        operator.convexhull = [];
        operator.bounds = @bounds;
        operator.derivative = @(x)(1./(x.*(x.^2-1).^0.5));

        varargout{1} = [];
        varargout{2} = operator;
        varargout{3} = varargin{3};

    otherwise
        error('SDPVAR/ASEC called with CHAR argument?');
end

function [L,U] = bounds(xL,xU)
if xL<=-1 & xU >=1
    L = -inf;
    U = inf;
else
    L = asec(xU);
    U = asec(xL);
end