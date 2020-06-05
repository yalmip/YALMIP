function varargout = asech(varargin)

switch class(varargin{1})

    case 'sdpvar'
        varargout{1} = InstantiateElementWise(mfilename,varargin{:});

    case 'char'

        operator = CreateBasicOperator('positive','callback');        
        operator.monotonicity = @monotonicity;        
        operator.bounds = @bounds;        
        operator.derivative = @(x)(1./(x.*(x.^2-1).^0.5));

        varargout{1} = [];
        varargout{2} = operator;
        varargout{3} = varargin{3};

    otherwise
        error(['SDPVAR/' upper(mfilename) ' called with weird argument']);
end

function mono = monotonicity(xL,xU)
if xU <= 0
    mono = 'increasing';
elseif xL >= 0
    mono = 'decreasing';
else
    mono = 'none';
end

function [L, U] = bounds(xL,xU)
if xU <= 0
    L = asech(xL);
    U = asech(xU);
elseif xL >= 0
    U = asech(xL);
    L = asech(xU);
else
    L = min(asech(xL),asech(xU));
    U = inf;
end