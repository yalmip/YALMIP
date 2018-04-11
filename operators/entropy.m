function varargout = entropy(varargin)
%ENTROPY
%
% y = ENTROPY(x)
%
% Computes/declares concave entropy -sum(x.*log(x))
%
% Implemented as evalutation based nonlinear operator. Hence, the concavity
% of this function is exploited to perform convexity analysis and rigorous
% modelling.
%
% See also CROSSENTROPY, KULLBACKLEIBLER.

switch class(varargin{1})

    case 'double'
        x = varargin{1};
        % Safe version with defined negative values (helps fmincon when
        % outside feasible region)
        if any(x<=0)
            z = abs(x);z(z==0)=1;
            y = z.*log(z);
            y(x==0) = 0;
            y(x<0) = 3*(x(x<0)-1).^2-3;
            varargout{1} = -sum(y);
        else
            varargout{1} = -sum(x.*log(x));
        end

    case {'sdpvar','ndsdpvar'}
 
        varargin{1} = reshape(varargin{1},[],1);
        varargout{1} = yalmip('define',mfilename,varargin{1});        

    case 'char'

        X = varargin{3};
        F = (X >= 0);

        operator = struct('convexity','concave','monotonicity','none','definiteness','none','model','callback');
        operator.range = [-inf exp(-1)*length(X)];
        operator.domain = [0 inf];
        operator.bounds = @bounds;
        operator.convexhull = @convexhull;
        operator.derivative = @derivative;
        varargout{1} = F;
        varargout{2} = operator;
        varargout{3} = X;

    otherwise
        error('SDPVAR/ENTROPY called with CHAR argument?');
end

function df = derivative(x)
x(x<=0)=eps;
df = (-1-log(x));

function [L, U] = bounds(xL,xU)

t = find(xL==0);
u = find(xL<0);
xL(t)=1;
LU = [-xL.*log(xL) -xU.*log(xU)];
LU(t,1)=0;
xL(t) = 0;
L = min(LU,[],2);
U = max(LU,[],2);
U((xL < exp(-1)) & (xU > exp(-1))) = exp(-1);
L = sum(L);
U = sum(U);


function [Ax, Ay, b] = convexhull(xL,xU)

if length(xL)==1
    xM = (xU+xL)/2;
    f1 = entropy(xL);
    f2 = entropy(xM);
    f3 = entropy(xU);        
    df1 = derivative(xL);
    df2 = derivative(xM);
    df3 = derivative(xU);
    [Ax,Ay,b] = convexhullConcave(xL,xM,xU,f1,f2,f3,df1,df2,df3);

elseif length(xL)==2
    x1 = [xL(1);xL(2)];
    x2 = [xU(1);xL(2)];
    x3 = [xL(1);xU(2)];
    x4 = [xU(1);xU(2)];
    x5 = (xL+xU)/2;
    
    f1 = entropy(x1);
    f2 = entropy(x2);
    f3 = entropy(x3);
    f4 = entropy(x4);
    f5 = entropy(x5);
    
    df1 = derivative(x1);
    df2 = derivative(x2);
    df3 = derivative(x3);
    df4 = derivative(x4);
    df5 = derivative(x5);
    
    [Ax,Ay,b] = convexhullConcave2D(x1,f1,df1,x2,f2,df2,x3,f3,df3,x4,f4,df4,x5,f5,df5);
else
    Ax = [];
    Ay = [];
    b = [];
end