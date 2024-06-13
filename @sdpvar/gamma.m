function varargout = gamma(varargin)

switch class(varargin{1})

    case 'sdpvar'
        varargout{1} = InstantiateElementWise(mfilename,varargin{:});

    case 'char'

        operator = CreateBasicOperator('convex','positive','callback');                                  
        operator.derivative =@(x)psi(0,x).*gamma(x);  
        operator.stationary = [1.46163214496836234126 8.856031944108888e-01];
        operator.range = [gamma(1.46163214496836234126) inf];
        operator.domain = [0 inf];
        
        varargout{1} = [];
        varargout{2} = operator;
        varargout{3} = varargin{3};

    otherwise
        error(['SDPVAR/' upper(mfilename) ' called with weird argument']);
end