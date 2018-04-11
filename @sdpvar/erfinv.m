function varargout = erfinv(varargin)
%ERFINV (overloaded)

switch class(varargin{1})

    case 'double'
        error('Overloaded SDPVAR/ERFINV CALLED WITH DOUBLE. Report error')

    case 'sdpvar'
        varargout{1} = InstantiateElementWise(mfilename,varargin{:});

    case 'char'
       
        X = varargin{3};
        F = (-1+1e-9 <= X <= 1-1e-9);

        operator = struct('convexity','none','monotonicity','increasing','definiteness','none','model','callback');
        operator.bounds = @bounds;

        varargout{1} = F;
        varargout{2} = operator;
        varargout{3} = X;

    otherwise
        error('SDPVAR/ERF called with CHAR argument?');
end

function [L,U] = bounds(xL,xU)
if xL<=-1
    L = -inf
else
    L = erfinv(xL);
end
if xU>=1
    U = inf;
else
    U = erfinv(xU);
end