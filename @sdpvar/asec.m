function varargout = asec(varargin)

switch class(varargin{1})

    case 'sdpvar'
        varargout{1} = InstantiateElementWise(mfilename,varargin{:});

    case 'char'

        operator = CreateBasicOperator('callback');        
        operator.bounds = @bounds;
        operator.derivative = @(x)(1./(x.*(x.^2-1).^0.5));

        varargout{1} = [];
        varargout{2} = operator;
        varargout{3} = varargin{3};

    otherwise
        error(['SDPVAR/' upper(mfilename) ' called with weird argument']);
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