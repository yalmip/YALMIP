function varargout = gamma(varargin)

switch class(varargin{1})

    case 'sdpvar'
        varargout{1} = InstantiateElementWise(mfilename,varargin{:});

    case 'char'

        operator = CreateBasicOperator('convex','positive','callback')               
        operator.monotonicity = @monotonicity;       
        operator.range = [gamma(1.46163214496836234126) inf];
        operator.derivative =@(x)psi(0,x).*gamma(x);  
        operator.stationary = [8.856031944108888e-01 1.46163214496836234126];
        
        varargout{1} = [];
        varargout{2} = operator;
        varargout{3} = varargin{3};

    otherwise
        error(['SDPVAR/' upper(mfilename) ' called with weird argument']);
end

function mono = monotonicity(xL,xU)
m = 1.46163214496836234126;
if xL < m && xU > m
    mono = 'none';
elseif xL >= m
    mono = 'increasing';
else
    mono = 'decreasing';
end