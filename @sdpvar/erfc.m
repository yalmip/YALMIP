function varargout = erfc(varargin)
%ERFC (overloaded)

switch class(varargin{1})

    case 'double'
        error('Overloaded SDPVAR/ERFC CALLED WITH DOUBLE. Report error')

    case 'sdpvar'
        varargout{1} = InstantiateElementWise(mfilename,varargin{:});

    case 'char'

        operator = struct('convexity',@convexity,'monotonicity','decreasing','definiteness','positive','model','callback');
        operator.bounds = @bounds;
        operator.range = [0 2];
        operator.derivative =@(x)-exp(-x.^2)*2/sqrt(pi);
        operator.inverse = @(x)(erfcinv(x));
         
        varargout{1} = [];
        varargout{2} = operator;
        varargout{3} = varargin{3};

    otherwise
        error('SDPVAR/ERF called with CHAR argument?');
end

function [L,U] = bounds(xL,xU)
L = erfc(xU);
U = erfc(xL);

function vexity = convexity(xL,xU)
if xL >= 0
    vexity = 'convex';
elseif xU <= 0
    vexity = 'concave';
else
    vexity = 'none';
end