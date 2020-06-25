function varargout = acot(varargin)

switch class(varargin{1})

    case 'sdpvar'
        varargout{1} = InstantiateElementWise(mfilename,varargin{:});

    case 'char'

        operator = CreateBasicOperator('callback');  
        operator.monotonicity = @monotonicity;   
        operator.bounds = @bounds;
        operator.derivative = @(x)(-(1 + x.^2).^-1);
        operator.range = [-pi/2 pi/2];
        operator.inflection = [0 1];

        varargout{1} = [];
        varargout{2} = operator;
        varargout{3} = varargin{3};

    otherwise
        error(['SDPVAR/' upper(mfilename) ' called with weird argument']);
end

function [L,U] = bounds(xL,xU)
if xL<=0 & xU >=0
    L = -pi/2;
    U = pi/2;
elseif xL>=0
    L = acot(xU);
    U = acot(xL);
else
    L = acot(xU);
    U = acot(xL);    
end

function mono = monotonicity(xL,xU)
if xL >= 0 || xU <= 0
    mono = 'decreasing';
else
    mono = 'none';
end