function varargout = erfcinv(varargin)

switch class(varargin{1})
    
    case 'sdpvar'
        varargout{1} = InstantiateElementWise(mfilename,varargin{:});
        
    case 'char'
        
        operator = CreateBasicOperator('decreasing','callback');
        operator.inverse = @(x)(erfc(x));
        operator.derivative = @(x)-1./(exp(-((erfcinv(x))).^2)*2/sqrt(pi));
        operator.inflection = [1 -1];
        operator.domain = [0 2];
        
        varargout{1} = [];
        varargout{2} = operator;
        varargout{3} = varargin{3};
        
    otherwise
        error(['SDPVAR/' upper(mfilename) ' called with weird argument']);
end