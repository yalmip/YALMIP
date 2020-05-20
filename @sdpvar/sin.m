function varargout = sin(varargin)
%SIN (overloaded)

switch class(varargin{1})

    case 'double'
        error('Overloaded SDPVAR/SIN CALLED WITH DOUBLE. Report error')

    case 'sdpvar'
        varargout{1} = InstantiateElementWiseUnitary(mfilename,varargin{:});
        %varargout{1} = InstantiateElementWise(mfilename,varargin{:});

    case 'char'

        % General operator
        operator = struct('monotonicity','none',...
            'definiteness','none',...
            'model','callback');

        operator.convexity  = @convexity;
        operator.bounds     = @bounds;     
        operator.derivative = @(x)(cos(x));
        operator.range = [-1 1];
        operator.domain = [-inf inf];     
  
        varargout{1} = [];
        varargout{2} = operator;
        varargout{3} = varargin{3};

    otherwise
        error('SDPVAR/SIN called with CHAR argument?');
end

function vexity = convexity(xL,xU)

if sin(xL)>=0 & sin(xU)>=0 & xU-xL<pi
    vexity = 'concave';    
elseif sin(xL)<=0 & sin(xU)<=0 & xU-xL<pi
    vexity = 'convex';   
else
    vexity = 'none';
end

function [L,U] = bounds(xL,xU)
if xU-xL >= 2*pi
    L = -1;
    U = 1;
else
    n = floor(( (xL + xU)/2/(2*pi)));
    xL = xL - n*2*pi;
    xU = xU - n*2*pi;
    yL = sin(xL);
    yU = sin(xU);
    L = min([yL yU]);
    U = max([yL yU]);
    if (xL<pi/2 & xU>pi/2) |  (xL<-3*pi/2 & xU>-3*pi/2)
        U = 1;
    end
    if (xL < 3*pi/2 & xU > 3*pi/2) | (xL < -pi/2 & xU > -pi/2)
        L = -1;
    end
end