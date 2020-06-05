function varargout = erf(varargin)

switch class(varargin{1})

    case 'sdpvar'
        varargout{1} = InstantiateElementWise(mfilename,varargin{:});

    case 'char'

        operator = CreateBasicOperator('increasing','callback');                      
        operator.derivative =@(x)exp(-x.^2)*2/sqrt(pi);
        operator.inverse = @(x)(erfinv(x));
        operator.inflection = [0 -1];
        operator.range = [-1 1];
              
        varargout{1} = [];
        varargout{2} = operator;
        varargout{3} = varargin{3};

    otherwise
        error(['SDPVAR/' upper(mfilename) ' called with weird argument']);
end