function varargout = xexpintinv(varargin)
%XEXPINTINV EXPINT(1/Z)/Z

switch class(varargin{1})

    case 'double'
        z = varargin{1};
        varargout{1} = (1./z).*expint(1./z);
        
    case 'sdpvar'
        varargout{1} = InstantiateElementWise(mfilename,varargin{:});

    case 'char'
        
        operator = CreateBasicOperator('positive','callback');
        operator.derivative = @derivative;
        operator.range = [0 1];
        operator.domain = [0 inf];
        
        varargout{1} = [];
        varargout{2} = operator;
        varargout{3} = varargin{3};

    otherwise
        error([upper(mfilename) ' called with weird argument']);
end

function d = derivative(z);
d = (-1./z.^2).*expint(1./z)+(1./z).*(-z.*exp(-1./z)).*(-1./z.^2);
