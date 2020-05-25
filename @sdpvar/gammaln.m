function varargout = gammaln(varargin)
%GAMMALN (overloaded)

switch class(varargin{1})

    case 'sdpvar'
        varargout{1} = InstantiateElementWise(mfilename,varargin{:});

    case 'char'

        operator = CreateBasicOperator('convex','callback');       
        operator.monotonicity = @monotonicity;        
        operator.domain = [0 inf];  
        operator.range = [gammaln(1.46163214496836234126) inf];
        operator.derivative = @(x)psi(0,x);       
        operator.stationary = [-1.214862905358496e-01 1.46163214496836234126];
        
        varargout{1} = [];
        varargout{2} = operator;
        varargout{3} = varargin{3};

    otherwise
        error('SDPVAR/GAMMALN called with CHAR argument?');
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