function varargout = plog(varargin)
%PLOG
%
% y = PLOG(x)
%
% Computes concave perspective log, x(1)*log(x(2)/x(1)) on x>0
%
% Implemented as evalutation based nonlinear operator. Hence, the concavity
% of this function is exploited to perform convexity analysis and rigorous
% modelling.

switch class(varargin{1})
    
    case 'double'
        
        if ~isequal(prod(size(varargin{1})),2)
            error('PLOG only defined for 2x1 arguments');
        end
        x = varargin{1};
        % Safe version with defined negative values (helps fmincon when
        % outside feasible region)

        if isequal(x(1),[0])
            varargout{1} = 0;
        else
            varargout{1} = x(1)*log(x(2)/x(1));
        end

    case 'sdpvar'

        if ~isequal(prod(size(varargin{1})),2)
            error('PLOG only defined for 2x1 arguments');
        else
            varargout{1} = yalmip('define',mfilename,varargin{1});
        end

    case 'char'

        X = varargin{3};
      
        operator = struct('convexity','concave','monotonicity','none','definiteness','none','model','callback');
        operator.range = [-inf inf];
        operator.domain = [0 inf];
        operator.bounds = @bounds;
        operator.convexhull = @convexhull;
        operator.derivative = @derivative;

        varargout{1} = [];
        varargout{2} = operator;
        varargout{3} = X;

    otherwise
        error('SDPVAR/PLOG called with CHAR argument?');
end

function dp = derivative(x)
z = x(2)/x(1);
dp = [log(z)-1;1./z];

function [L,U] = bounds(xL,xU)
xU(isinf(xU)) = 1e12;
x1 = xL(1)*log(xL(2)/xL(1));
x2 = xU(1)*log(xU(2)/xU(1));
x3 = xL(1)*log(xU(2)/xL(1));
x4 = xU(1)*log(xL(2)/xU(1));
L = min([x1 x2 x3 x4]);

% Stationary in x1 when x1 = x2*exp(-1)
% Increasing in x2, so max at border
p1 = [exp(-1)*xU(2);xU(2)];
x5 = p1(1)*log(p1(2)/p1(1));
if x5(1) >= xL(1) && x5(1) <= xU(1)    
    U = max([x1 x2 x3 x4 x5]);
else
    U = max([x1 x2 x3 x4]);
end

function [Ax,Ay,b] = convexhull(xL,xU)

x1 = [xL(1);xL(2)];
x2 = [xU(1);xL(2)];
x3 = [xL(1);xU(2)];
x4 = [xU(1);xU(2)];
x5 = (xL+xU)/2;

f1 = plog(x1);
f2 = plog(x2);
f3 = plog(x3);
f4 = plog(x4);
f5 = plog(x5);

df1 = derivative(x1);
df2 = derivative(x2);
df3 = derivative(x3);
df4 = derivative(x4);
df5 = derivative(x5);

[Ax,Ay,b] = convexhullConcave2D(x1,f1,df1,x2,f2,df2,x3,f3,df3,x4,f4,df4,x5,f5,df5);
