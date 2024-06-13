function varargout = erfc(varargin)

switch class(varargin{1})

    case 'sdpvar'
        varargout{1} = InstantiateElementWise(mfilename,varargin{:});

    case 'char'

        operator = CreateBasicOperator('decreasing','positive','callback');                      
        operator.derivative =@(x)-exp(-x.^2)*2/sqrt(pi);
        operator.inverse = @(x)(erfcinv(x));
        operator.inflection = [-inf -1 0 1];
        operator.range = [0 2];
         
        varargout{1} = [];
        varargout{2} = operator;
        varargout{3} = varargin{3};

    otherwise
        error(['SDPVAR/' upper(mfilename) ' called with weird argument']);
end