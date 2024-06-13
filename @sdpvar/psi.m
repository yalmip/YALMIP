function varargout = psi(varargin)

switch class(varargin{1})

    case 'sdpvar'
        if nargin > 1
            error('PSI currently only supported with 1 argument. Request feature if needed');
        end
        varargout{1} = InstantiateElementWise(mfilename,varargin{:});

    case 'char'

        operator = CreateBasicOperator('concave','increasing','callback');               
        operator.domain = [0 inf];  
        operator.range = [0 700];
        operator.derivative = @(x)psi(1,x);       
                
        varargout{1} = [];
        varargout{2} = operator;
        varargout{3} = varargin{3};

    otherwise
        error(['SDPVAR/' upper(mfilename) ' called with weird argument']);
end