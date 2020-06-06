function varargout = sech(varargin)

switch class(varargin{1})

    case 'sdpvar'
        varargout{1} = InstantiateElementWise(mfilename,varargin{:});

    case 'char'

        operator = CreateBasicOperator('callback');
        operator.monotonicity = @monotonicity;                
        operator.derivative = @(x)(-tanh(x).*sech(x)); 
        operator.stationary = [0 1];
        operator.inflection = [-0.881493604 -1 0.881493604 1];
        operator.range = [0 1];
        
        varargout{1} = [];
        varargout{2} = operator;
        varargout{3} = varargin{3};

    otherwise
        error(['SDPVAR/' upper(mfilename) ' called with weird argument']);
end

function mono = monotonicity(xL,xU)
if xL >= 0  
    mono = 'decreasing';
elseif xU <= 0
    mono = 'increasing';
else
    mono = 'none';
end