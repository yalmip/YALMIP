function varargout = cdf_internal(varargin)

switch class(varargin{1})
    
    case 'double'
        varargout{1} = cdf(varargin{2},varargin{1},varargin{3:end});
        
    case 'char'
        
        operator = struct('convexity','none','monotonicity','increasing','definiteness','positive','model','callback');
        operator.bounds = @bounds;
        
        varargout{1} = [];
        varargout{2} = operator;
        varargout{3} = varargin{3};
        
    otherwise
        error('SDPVAR/ERF called with CHAR argument?');
end

function [L,U] = bounds(xL,xU,varargin)
L = cdf(varargin{1},xL,varargin{2:end});
U = cdf(varargin{1},xU,varargin{2:end});