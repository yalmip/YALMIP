function varargout = cdf_internal(varargin)

switch class(varargin{1})
    
    case 'double'
        varargout{1} = cdf(varargin{2},varargin{1},varargin{3:end});
        
    case 'char'
        
        operator = struct('convexity','none','monotonicity','increasing','definiteness','positive','model','callback');
        operator.bounds = @bounds;
        operator.range = [0 1];
        operator.derivative = @(x)derivative(x,varargin{4:end});
        operator.inverse = @(x)icdf(varargin{4},x,varargin{5:end});  
        
        varargout{1} = [];
        varargout{2} = operator;
        varargout{3} = varargin{3};
        
    otherwise
        error('SDPVAR/CDF called with CHAR argument?');
end

function df = derivative(x,varargin)
df = pdf(varargin{1},x,varargin{2:end});

function [L,U] = bounds(xL,xU,varargin)
L = cdf(varargin{1},xL,varargin{2:end});
U = cdf(varargin{1},xU,varargin{2:end});