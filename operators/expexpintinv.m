function varargout = expexpintinv(varargin)

switch class(varargin{1})

    case 'double'
        z = varargin{1};        
        varargout{1} = exp(1./z).*expint(1./z);
        
    case 'sdpvar'
        varargout{1} = InstantiateElementWise(mfilename,varargin{:});

    case 'char'
        
        operator = CreateBasicOperator('concave','increasing','positive','callback');
        operator.derivative = @derivative;
        operator.range = [0 inf];
        operator.domain = [0 inf];
        
        varargout{1} = [];
        varargout{2} = operator;
        varargout{3} = varargin{3};

    otherwise
        error([upper(mfilename) ' called with weird argument']);
end

function d = derivative(z);
d = (-1./z.^2).*exp(1./z).*expint(1./z)+exp(1./z).*(-z.*exp(-1./z)).*(-1./z.^2);