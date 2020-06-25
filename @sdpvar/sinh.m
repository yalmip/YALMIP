function varargout = sinh(varargin)

switch class(varargin{1})

    case 'sdpvar'
        varargout{1} = InstantiateElementWise(mfilename,varargin{:});
 
    case 'char'

        operator = CreateBasicOperator('increasing','odd','callback');        
        operator.derivative = @(x)(cosh(x));
        operator.inflection = [0 1];

        varargout{1} = [];
        varargout{2} = operator;
        varargout{3} = varargin{3};

    otherwise
        error(['SDPVAR/' upper(mfilename) ' called with weird argument']);
end