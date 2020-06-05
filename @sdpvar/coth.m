function varargout = coth(varargin)

switch class(varargin{1})

    case 'sdpvar'
        varargout{1} = InstantiateElementWise(mfilename,varargin{:});

    case 'char'

        operator = CreateBasicOperator('callback');
        operator.monotonicity = @monotonicity;                
        operator.convexity = @convexity;                
        operator.derivative = @(x)(-(1-coth(x)).^2);        
               
        varargout{1} = [];
        varargout{2} = operator;
        varargout{3} = varargin{3};

    otherwise
        error(['SDPVAR/' upper(mfilename) ' called with weird argument']);
end

function mono = monotonicity(xL,xU)
if xL >= 0 || xU <= 0
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