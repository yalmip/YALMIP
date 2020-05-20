function varargout = atan(varargin)
%ATAN (overloaded)

switch class(varargin{1})

    case 'sdpvar'
        varargout{1} = InstantiateElementWise(mfilename,varargin{:});

    case 'char'

        operator = struct('convexity',@convexity,'monotonicity','increasing','definiteness','none','model','callback');        
        operator.bounds = @bounds;
        operator.derivative = @(x)((1+x.^2).^-1);
        operator.inverse = @(x)(tan(x));
        operator.range = [-pi/2 pi/2];
        
        varargout{1} = [];
        varargout{2} = operator;
        varargout{3} = varargin{3};

    otherwise
        error('SDPVAR/ATAN called with CHAR argument?');
end

function [L,U] = bounds(xL,xU)
L = atan(xL);
U = atan(xU);

function vexity = convexity(xL,xU)
if xL >= 0  
    vexity = 'concave';
elseif xU <= 0
    vexity = 'convex';
else
    vexity = 'none';
end