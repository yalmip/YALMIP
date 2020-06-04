function varargout = expint(varargin)
%EXPINT (overloaded)

switch class(varargin{1})

    case 'sdpvar'
        varargout{1} = InstantiateElementWise(mfilename,varargin{:});

    case 'char'
        
        operator = CreateBasicOperator('convex','decreasing','positive','callback');        
        operator.derivative = @(x)(-exp(-x)./x);
        operator.domain = [0 inf];
        varargout{1} = [];
        varargout{2} = operator;
        varargout{3} = varargin{3};

    otherwise
        error('SDPVAR/EXPINT called with CHAR argument?');
end
