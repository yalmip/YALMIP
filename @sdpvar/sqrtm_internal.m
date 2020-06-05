function varargout = sqrtm_internal(varargin)

switch class(varargin{1})

    case {'double', 'gem', 'sgem'}
        varargout{1} = sqrt(varargin{1});

    case 'sdpvar'
          varargout{1} = InstantiateElementWise(mfilename,varargin{:});

    case 'char'
                        
        operator = CreateBasicOperator('concave','increasing','positive','callback');
        operator.derivative = @(x) 1./(eps + 2*abs(x).^0.5);
        operator.inverse = @(x)x.^2;
        operator.domain = [0 inf];

        varargout{1} = [];
        varargout{2} = operator;
        varargout{3} = varargin{3};

    otherwise
        error(['SDPVAR/' upper(mfilename) ' called with weird argument']);
end