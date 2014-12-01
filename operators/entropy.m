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
% See also crossentropy.

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

    case 'sdpvar'

        if min(size(varargin{1}))>1
            error('ENTROPY only defined for vector arguments');
        else
            varargout{1} = yalmip('define',mfilename,varargin{1});
        end

    case 'char'

        X = varargin{3};
        F = (X >= 0);

        operator = struct('convexity','concave','monotonicity','none','definiteness','none','model','callback');
        operator.range = [-inf exp(-1)*length(X)];
        operator.domain = [0 inf];
        operator.bounds = @bounds;
        operator.convexhull = @convexhull;
        operator.derivative = @(x) (-1-log(x));

        varargout{1} = F;
        varargout{2} = operator;
        varargout{3} = X;

    otherwise
        error('SDPVAR/LOG called with CHAR argument?');
end


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
Ax = [];
Ay = [];
b = [];
% % Loop thorough all variables, and compute a convex hull each term xlogx
% Hmm, how do I merge without expensive projection
% for i = 1:length(xL)
%     if xL(i) <= 0
%         fL = inf;
%         dfL = -inf;
%     else
%         fL = -xL(i).*log(xL(i));
%         dfL = -log(xL(i)) - 1;
%     end
%     fU = -xU(i).*log(xU(i));
%     dfU = -log(xU(i)) - 1;
%     [Ax,Ay,b] = convexhullConcave(xL(i),xU(i),fL,fU,dfL,dfU);
% end