function varargout = erf(varargin)
%ERF (overloaded)

switch class(varargin{1})

    case 'double'
        error('Overloaded SDPVAR/ERF CALLED WITH DOUBLE. Report error')

    case 'sdpvar'
        varargout{1} = InstantiateElementWise(mfilename,varargin{:});

    case 'char'

        operator = struct('convexity',@convexity,'monotonicity','increasing','definiteness','none','model','callback');
        operator.bounds = @bounds;
        operator.range = [-1 1];
        operator.derivative =@(x)exp(-x.^2)*2/sqrt(pi);
        operator.inverse = @(x)(erfinv(x));
              
        varargout{1} = [];
        varargout{2} = operator;
        varargout{3} = varargin{3};

    otherwise
        error('SDPVAR/ERF called with CHAR argument?');
end

function [L,U] = bounds(xL,xU)
L = erf(xL);
U = erf(xU);

function vexity = convexity(xL,xU)
if xL >= 0  
    vexity = 'concave';
elseif xU <= 0
    vexity = 'convex';
else
    vexity = 'none';
end

