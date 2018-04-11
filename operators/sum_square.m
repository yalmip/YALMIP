function varargout = sum_square(varargin)

switch class(varargin{1})

    case 'double'
        x = varargin{1};
        varargout{1} = x(:)'*x(:);
        
    case 'sdpvar'         
        varargout{1} = yalmip('define',mfilename,varargin{1});        

    case 'char'

        operator = struct('convexity','convex','monotonicity','none','definiteness','positive','model','callback');
        operator.range = [0 inf];
        operator.domain = [-inf inf];                
        operator.derivative = @(x)(2*x);
        operator.bounds = @bounds;
        varargout{1} = [];
        varargout{2} = operator;
        varargout{3} = varargin{3};

    otherwise
        error('SDPVAR/SUM_SQUARE called with strange argument.');
end

function [L,U] = bounds(xL,xU);
pass0 = xU > 0 & xL < 0;
low   = min([xL xU],[],2).^2;
high  = max([xL xU],[],2).^2;
low(pass0) = 0;
L = sum(low);
U = sum(high);
