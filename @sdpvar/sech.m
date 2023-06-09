function varargout = sech(varargin)

switch class(varargin{1})

    case 'sdpvar'
        varargout{1} = InstantiateElementWise(mfilename,varargin{:});

    case 'char'

        operator = CreateBasicOperator('bell-shape','callback');        
        operator.derivative = @(x)(-tanh(x).*sech(x)); 
        operator.stationary = [0 1];
        operator.inflection = [-inf 1 -0.881493604 -1 0.881493604 1];
        operator.range = [0 1];
        
        varargout{1} = [];
        varargout{2} = operator;
        varargout{3} = varargin{3};

    otherwise
        error(['SDPVAR/' upper(mfilename) ' called with weird argument']);
end