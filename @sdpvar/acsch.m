function varargout = acsch(varargin)

switch class(varargin{1})

    case 'sdpvar'
        varargout{1} = InstantiateElementWise(mfilename,varargin{:});

    case 'char'

        operator = CreateBasicOperator('callback');
        operator.monotinicity = @monotinicity;
        operator.derivative = @(x)(-1./(abs(x).*sqrt(1 + x.^2)));
        operator.singularity = 0;
        operator.inflection = [0 1];

        varargout{1} = [];
        varargout{2} = operator;
        varargout{3} = varargin{3};

    otherwise
        error(['SDPVAR/' upper(mfilename) ' called with weird argument']);
end

function mono = monotinicity(xL,xU)
if xU <= 0
    mono = 'decreasing';
elseif xL >= 0
    mono = 'decreasing';
else
    mono = 'none';
end

function vexity = convexity(xL,xU)
if xL >= 0  
    vexity = 'convex';
elseif xU <= 0
    vexity = 'concave';
else
    vexity = 'none';
end

