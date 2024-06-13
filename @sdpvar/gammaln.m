function varargout = gammaln(varargin)

switch class(varargin{1})

    case 'sdpvar'
        varargout{1} = InstantiateElementWise(mfilename,varargin{:});

    case 'char'

        operator = CreateBasicOperator('convex','callback');                     
        operator.domain = [0 inf];  
        operator.range = [gammaln(1.46163214496836234126) inf];
        operator.derivative = @(x)psi(0,x);       
        operator.stationary = [1.46163214496836234126 -1.214862905358496e-01];
        
        varargout{1} = [];
        varargout{2} = operator;
        varargout{3} = varargin{3};

    otherwise
        error(['SDPVAR/' upper(mfilename) ' called with weird argument']);
end