function varargout = pexp(varargin)
%PEXP
%
% y = PEXP(x)
%
% Defines perspective exp, x(1)*exp(x(2)/x(1)) on x(1)>0
% 
% Alternatively
%
% y = PEXP(x,y) to define x*exp(y/x) on x>0
%
% Implemented as either evalutation based-nonlinear operator, or
% represented using exponential cones depending on solver. Hence, the
% convexity of this function is exploited to perform convexity analysis and
% rigorous modelling. 
%
% See also ENTROPY, LOGSUMEXP, CROSSENTROPY, KULLBACKLEIBLER, EXPCONE

switch class(varargin{1})
    
    case 'double'
        
        if nargin == 2
            varargin{1} = [varargin{1};varargin{2}];
        end
        if nargin == 1 && ~isequal(prod(size(varargin{1})),2)
            error('PEXP only defined for 2x1 arguments');
        end
        x = varargin{1};
        
        varargout{1} = x(1)*exp(x(2)/x(1));


    case 'sdpvar'
        
        if nargin == 2
            varargin{1} = [varargin{1};varargin{2}];
        end               
        if ~isequal(prod(size(varargin{1})),2)
            error('PEXP only defined for 2x1 arguments');
        else
            varargout{1} = yalmip('define',mfilename,varargin{1});
        end

    case 'char'
        
        operator = CreateBasicOperator('convex','positive','callback');
        operator.range = [0 inf];   
        operator.derivative = @derivative;
        operator.convexhull = @convexhull;
        operator.bounds = @bounds;
        
        varargout{1} = [varargin{3}(1) >= 0];
        varargout{2} = operator;
        varargout{3} = varargin{3};

    otherwise
        error([upper(mfilename) ' called with weird argument']);
end

function [L,U] = bounds(xL,xU)
x1 = [xL(1);xL(2)];
x2 = [xU(1);xL(2)];
x3 = [xL(1);xU(2)];
x4 = [xU(1);xU(2)];
L = min([pexp(x1) pexp(x2) pexp(x3) pexp(x4)]);
U = max([pexp(x1) pexp(x2) pexp(x3) pexp(x4)]);

function dp = derivative(x)
z = x(2)/x(1);
dp = [exp(z)-z*exp(z);exp(z)];

function [Ax,Ay,b,K] = convexhull(xL,xU)
x1 = [xL(1);xL(2)];
x2 = [xU(1);xL(2)];
x3 = [xL(1);xU(2)];
x4 = [xU(1);xU(2)];
x5 = (xL+xU)/2;
f1 = pexp(x1);
f2 = pexp(x2);
f3 = pexp(x3);
f4 = pexp(x4);
f5 = pexp(x5);
df1 = derivative(x1);
df2 = derivative(x2);
df3 = derivative(x3);
df4 = derivative(x4);
df5 = derivative(x5);
[Ax,Ay,b,K] = convexhullConvex2D(x1,f1,df1,x2,f2,df2,x3,f3,df3,x4,f4,df4,x5,f5,df5);