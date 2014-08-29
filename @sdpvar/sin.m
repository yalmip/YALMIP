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
        operator = struct('convexity','none',...
            'monotonicity','none',...
            'definiteness','none',...
            'model','callback');

        operator.bounds     = @bounds;
        operator.convexhull = @convexhull;
        operator.derivative = @(x)(cos(x));
        operator.range = [-1 1];
        operator.domain = [-inf inf];     
  
        varargout{1} = [];
        varargout{2} = operator;
        varargout{3} = varargin{3};

    otherwise
        error('SDPVAR/SIN called with CHAR argument?');
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

function [Ax, Ay, b] = convexhull(xL,xU)
if sin(xL)>=0 & sin(xU)>=0 & xU-xL<pi
    xM = (xL+xU)/2;
    fL = sin(xL);
    fM = sin(xM);
    fU = sin(xU);
    dfL = cos(xL);
    dfM = cos(xM);
    dfU = cos(xU);
    [Ax,Ay,b] = convexhullConcave(xL,xM,xU,fL,fM,fU,dfL,dfM,dfU);
elseif sin(xL)<=0 & sin(xU)<=0 & xU-xL<pi
    fL = sin(xL);
    fU = sin(xU);
    dfL = cos(xL);
    dfU = cos(xU);
    [Ax,Ay,b] = convexhullConvex(xL,xU,fL,fU,dfL,dfU);
else
    [Ax,Ay,b] = convexhullGeneral(xL,xU,@sin);
end




