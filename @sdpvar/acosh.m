function varargout = acosh(varargin)

switch class(varargin{1})

    case 'sdpvar'
        varargout{1} = InstantiateElementWise(mfilename,varargin{:});

    case 'char'

        operator = CreateBasicOperator('positive','callback');   
        operator.convexity = @convexity;
        operator.bounds = @bounds;

        varargout{1} = [];
        varargout{2} = operator;
        varargout{3} = varargin{3};

    otherwise
        error(['SDPVAR/' upper(mfilename) ' called with weird argument']);
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

function vexity = convexity(xL,xU)
if xL >= 1
    vexity = 'concave';
elseif xU <= -1
    vexity = 'concave';
else
    vexity = 'none';
end
