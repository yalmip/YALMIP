function varargout = icdf_internal(varargin)

switch class(varargin{1})
    
    case 'double'
        varargout{1} = icdf(varargin{2},varargin{1},varargin{3:end});
        
    case 'char'
        
        operator = struct('convexity','none','monotonicity','increasing','definiteness','none','model','callback');
        operator.bounds = @bounds;
        operator.domain = [1e-16 1-1e-16];
        operator.derivative = @(x)derivative(x,varargin{4:end});
        
        varargout{1} = [1-1e-8 >= varargin{3} >= 1e-8];
        varargout{2} = operator;
        varargout{3} = varargin{3};
        
    otherwise
        error('SDPVAR/ICDF called with CHAR argument?');
end

function df = derivative(x,varargin)

PsiInv = icdf(varargin{1},x,varargin{2:end});
Phi = pdf(varargin{1},PsiInv,varargin{2:end});
df = 1/Phi;

function [L,U] = bounds(xL,xU,varargin)
L = icdf(varargin{1},xL,varargin{2:end});
U = icdf(varargin{1},xU,varargin{2:end});