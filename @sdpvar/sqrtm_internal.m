function varargout = sqrtm_internal(varargin)
%SQRTM (overloaded)

switch class(varargin{1})

    case {'double', 'gem', 'sgem'}
        varargout{1} = sqrt(varargin{1});

    case 'sdpvar'
          varargout{1} = InstantiateElementWise(mfilename,varargin{:});

    case 'char'

        X = varargin{3};
        F = (X >= 0);
        
        operator = CreateBasicOperator('concave','increasing','positive','callback');
        operator.derivative = @(x) 1./(eps + 2*abs(x).^0.5);
        operator.inverse = @(x)x.^2;

        varargout{1} = F;
        varargout{2} = operator;
        varargout{3} = X;

    otherwise
        error(['SDPVAR/' upper(mfilename) ' called with weird argument']);
end