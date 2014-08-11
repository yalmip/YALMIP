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

        varargout{1} = [varargin{3}.^2 >= 1]; % Disconnected domain
        varargout{2} = operator;
        varargout{3} = varargin{3};

    otherwise
        error('SDPVAR/ASEC called with CHAR argument?');
end

function [L,U] = bounds(xL,xU)
if xU <= -1 || xL >= 1
    L = asec(xL);
    U = asec(xU);
elseif xL < 0 & xU > 0
    L = 0;
    U = pi;
elseif xU < 0 || xL > 0
    L = real(asec(xL));
    U = real(asec(xU));
else
    L = 0;
    U = pi;
end