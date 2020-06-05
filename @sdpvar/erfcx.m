function varargout = erfcx(varargin)

switch class(varargin{1})

    case 'sdpvar'
        varargout{1} = InstantiateElementWise(mfilename,varargin{:});

    case 'char'

        operator = CreateBasicOperator('convex','decreasing','positive','callback');        
        operator.derivative = @(x)(2*x.*exp(x.^2).*erfc(x) + exp(x.^2).*(-exp(-x.^2)*2/sqrt(pi)));

        varargout{1} = [];
        varargout{2} = operator;
        varargout{3} = varargin{3};

    otherwise
        error(['SDPVAR/' upper(mfilename) ' called with weird argument']);
end