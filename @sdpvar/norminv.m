function varargout = norminv(varargin)
%NORMINV (overloaded)

switch class(varargin{1})

    case 'double'
        error('Overloaded SDPVAR/NORMINV CALLED WITH DOUBLE. Report error')

    case 'sdpvar'
        varargout{1} = InstantiateElementWise(mfilename,varargin{:});

    case 'char'
       
        X = varargin{3};
        F = (0 <= X <= 1);

        operator = struct('convexity','none','monotonicity','increasing','definiteness','none','model','callback');
        operator.bounds = @bounds;

        varargout{1} = F;
        varargout{2} = operator;
        varargout{3} = X;

    otherwise
        error('SDPVAR/NORMINV called with CHAR argument?');
end

function [L,U] = bounds(xL,xU)
if xL<=0
    L = -inf
else
    L = norminv(xL);
end
if xU>=1
    U = inf;
else
    U = norminv(xU);
end