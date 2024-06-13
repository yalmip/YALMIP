function varargout = ellipke(varargin)

switch class(varargin{1})
    
    case 'sdpvar'
        varargout{1} = InstantiateElementWise(mfilename,varargin{:});
        
    case 'char'
        
        operator = CreateBasicOperator('convex','increasing','positive','callback')
        operator.derivative = @(x)(airy(1,x));
        operator.range = [1.570796326794897e+00 inf];
        
        varargout{1} = [];
        varargout{2} = operator;
        varargout{3} = varargin{3};
        
    otherwise
        error(['SDPVAR/' upper(mfilename) ' called with weird argument']);
end