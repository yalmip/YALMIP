function varargout = asinh(varargin)

switch class(varargin{1})

    case 'sdpvar'
        varargout{1} = InstantiateElementWise(mfilename,varargin{:});

    case 'char'

        operator = CreateBasicOperator('increasing','callback');                      
        operator.derivative = @(x)((1 + x.^2).^-0.5);
        operator.inflection = [0 -1];
        operator.range = [-700 700];
            
        varargout{1} = [];
        varargout{2} = operator;
        varargout{3} = varargin{3};

    otherwise
        error(['SDPVAR/' upper(mfilename) ' called with weird argument']);
end