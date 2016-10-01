function varargout = norm_callback(varargin)
% Low-level alternative for norm operations with callback solver

switch class(varargin{1})

    case 'double'
        x = varargin{1};        
        varargout{1} = sqrt(x(:)'*x(:));
        
    case {'sdpvar','ndsdpvar'} 
        varargin{1} = reshape(varargin{1},[],1);
        varargout{1} = yalmip('define',mfilename,varargin{1});        

    case 'char'

        X = varargin{3};
        F = [];

        operator = struct('convexity','convex','monotonicity','none','definiteness','positive','model','callback');
        operator.range = [0 inf];
        operator.domain = [-inf inf];
        operator.bounds = @bounds;
        operator.derivative = @derivative;
        varargout{1} = F;
        varargout{2} = operator;
        varargout{3} = X;

    otherwise
        error('NORM_CALLBACK called with weird argument.');
end

function df = derivative(x)
x = x(:);
df = x./(norm(x));
df(isinf(df))=0;
df(isnan(df))=0;

function [L, U] = bounds(xL,xU)
pass0 = xU > 0 & xL < 0;
low   = min([xL xU],[],2).^2;
high  = max([xL xU],[],2).^2;
low(pass0) = 0;
L = sqrt(sum(low));
U = sqrt(sum(high));
