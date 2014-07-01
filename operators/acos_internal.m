function varargout = acos_internal(varargin)
%ACOS (overloaded)

switch class(varargin{1})

    case 'double'
        x = varargin{1};
        y = acos(x);
        y(x<-1) = pi;
        y(x>1)  = 0;
        varargout{1} = y;
        
    case 'char'

        operator = struct('convexity','none','monotonicity','decreasing','definiteness','none','model','callback');
        operator.convexhull = [];
        operator.bounds = @bounds;
        operator.derivative = @derivative;
        operator.range = [-pi/2 pi/2];
        operator.domain = [-1 1];

        varargout{1} = [];
        varargout{2} = operator;
        varargout{3} = varargin{3};

    otherwise
        error('SDPVAR/ACOS called with CHAR argument?');
end

function [L,U] = bounds(xL,xU)
L = real(acos(xU));
U = real(acos(xL));

function df = derivative(x)

df = (-(1 - x.^2).^-0.5);
df(x>1) = 0;
df(x<1) = 0;
