function varargout = acoth(varargin)

switch class(varargin{1})

    case 'sdpvar'
        varargout{1} = InstantiateElementWise(mfilename,varargin{:});

    case 'char'

        operator = CreateBasicOperator('decreasing','convex','positive','callback');        
        operator.derivative = @(x)(1./(1-x.^2));        
        operator.inverse = @(x)(coth(x));        
        operator.domain = [1 inf];
        
        varargout{1} = [];
        varargout{2} = operator;
        varargout{3} = varargin{3};

    otherwise
        error(['SDPVAR/' upper(mfilename) ' called with weird argument']);
end