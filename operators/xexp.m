function varargout = xexp(varargin)
%XEXP X*exp(X)

switch class(varargin{1})

    case 'double'
        z = varargin{1};
        varargout{1} = z.*exp(z);
        
    case 'sdpvar'
        varargout{1} = InstantiateElementWise(mfilename,varargin{:});

    case 'char'
        
        operator = CreateBasicOperator('callback','convex');
        operator.derivative = @derivative;
        operator.range = [-1*exp(-1) inf];
        operator.domain = [-inf inf];
        
        varargout{1} = [];
        varargout{2} = operator;
        varargout{3} = varargin{3};

    otherwise
        error([upper(mfilename) ' called with weird argument']);
end

function d = derivative(z)
d = (1+z).*exp(z);
