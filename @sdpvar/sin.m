function varargout = sin(varargin)
%SIN (overloaded)

switch class(varargin{1})

    case 'sdpvar'
        varargout{1} = InstantiateElementWiseUnitary(mfilename,varargin{:});

    case 'char'
        
        operator = CreateBasicOperator('callback');
        operator.definiteness  = @definiteness;
        operator.monotonicity  = @monotonicity;
        operator.convexity  = @convexity;
        operator.bounds     = @bounds;     
        operator.stationary = @stationary; 
        operator.inflection = @inflection;
        operator.periodic = pi; 
        operator.derivative = @(x)(cos(x));
        operator.range = [-1 1];
               
        varargout{1} = [];
        varargout{2} = operator;
        varargout{3} = varargin{3};

    otherwise
        error(['SDPVAR/' upper(mfilename) ' called with weird argument']);
end

function def = definiteness(xL,xU)
if xU-xL > pi
    def = 'none';
else
     n = floor(xU/(2*pi));
    xL = xL - n*2*pi;
    xU = xU - n*2*pi;
	if xL >= 0 && xU <=pi
        def = 'positive';
    elseif (xL >= -pi && xU <= 0) || (xL >= pi && xU <= 2*pi)
        def = 'negative';
    else
        def = 'none';
    end
end

function mono = monotonicity(xL,xU)
if xU-xL > pi
    mono = 'none';
else
    n = floor(xU/(2*pi));
    xL = xL - n*2*pi;
    xU = xU - n*2*pi;
	if (xL >= -pi/2 && xU <=pi/2) || (xL >= 1.5*pi && xU <=2*pi)
        mono = 'increasing';
    elseif (xL >= -1.5*pi && xU <= -pi/2) || (xL >= pi/2 && xU <= 1.5*pi)
        mono = 'decreasing';
    else
        mono = 'none';
    end
end
    
function vexity = convexity(xL,xU)
if xU-xL > pi
    vexity = 'none';
else
    n = floor(xU/(2*pi));
    xL = xL - n*2*pi;
    xU = xU - n*2*pi;
	if (xL >= 0 && xU <=pi) || (xL >= -2*pi && xU <=-pi)
        vexity = 'concave';
    elseif (xL >= -pi && xU <= 0) || (xL >= pi && xU <= 2*pi)
        vexity = 'convex';
    else
        vexity = 'none';
    end
end

function inflections = inflection(xL,xU)
r = floor(xU/(pi/2));
t = ceil(xL/(pi/2));
spots = [t:r];
spots = spots(rem(spots,2)==0);
dir = double(abs(rem(spots/2,2))==1);
dir(dir==0)=-1;
xS = spots*pi/2;
if isempty(xS)
    inflections = [];
else
    xS = [-inf xS];
    dir = [-dir(1) dir];
    inflections = reshape([xS;dir],1,[]);
end
function [xS,fS] = stationary(xL,xU)
r = floor(xU/(pi/2));
t = ceil(xL/(pi/2));
spots = [t:r];
spots = spots(rem(spots,2)~=0);
xS = spots*pi/2;
fS = sign(sin(xS));

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