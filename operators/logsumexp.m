function varargout = logsumexp(varargin)
%LOGSUMEXP
%
% y = LOGSUMEXP(x)
%
% Computes/declares log of sum of exponentials log(sum(exp(x)))
%
% Implemented as evalutation based nonlinear operator. Hence, the convexity
% of this function is exploited to perform convexity analysis and rigorous
% modelling.

switch class(varargin{1})

    case 'double'
        x = varargin{1};
        varargout{1} = log(sum(exp(x)));

    case 'sdpvar'

        if min(size(varargin{1}))>1
            x = varargin{1};
            y = [];
            for i = 1:size(x,2)
                y = [y yalmip('define',mfilename,x(:,i))];
            end
            varargout{1} = y;   
        elseif max(size(varargin{1}))==1
            varargout{1} = varargin{1};
        else
            varargout{1} = yalmip('define',mfilename,varargin{1});
        end

    case 'char'
              
        operator = CreateBasicOperator('convex','callback');        
        operator.bounds = @bounds;
        operator.derivative = @(x) exp(x)./(sum(exp(x)));
        operator.convexhull = @convexhull;
                
        varargout{1} = [];
        varargout{2} = operator;
        varargout{3} = varargin{3};

    otherwise
        error([upper(mfilename) ' called with weird argument']);
end


function [L, U] = bounds(xL,xU)

L = log(sum(exp(xL)));
U = log(sum(exp(xU)));

function [Ax, Ay, b] = convexhull(xL,xU)
% Ax*x + Ay*y < b
% Use lower cut from center point
c = (xL + xU)/2;
df = exp(c)./sum(exp(c));
Ax = df';
Ay = -1;
b = -(log(sum(exp(c)))-df'*c);

% Add simple Fourier-Motzkins projections of disaggregated model
% log(sum(q)) upper bounded by y <= k*sum(q)+c, 
%                              q =(exp x)), 
%                              qi <= alpha*xi + gammai,
%  (y-c)/k <= sum(q), replace sum with sum of upper bounds
qL = sum(exp(xL));
qU = sum(exp(xU));
q = sum(exp(xL + linspace(0,1,3).*(xU-xL)),1);
% Pick one point in the middle. However, middle is not good as it adds
% nothing if xU is large. Without much though, pick point where gradient is
% average of gradient at xL and XU
qM = 2*(qU.*qL)./(qU + qL);
cL = log(qL) - (1/qL)*qL;
kL = (1/qL);
cU = log(qU) - (1/qU)*qU;
kU = (1/qU);
cM = log(qM) - (1/qM)*qM;
kM = (1/qM);
c = log(q) - (1./q).*q;
k = 1./q;

gamma = exp(xL) - (xL./(xU-xL)).*(exp(xU)-exp(xL));
alpha = (exp(xU)-exp(xL))./(xU-xL);

Ax = [Ax;-(k.*alpha)'];
Ay = [Ay;ones(length(q),1)];
b = [b;(c+k*sum(gamma))'];

% Add a cut based on sum(x) <= prod x when x >= n^(1/(n-1))
n=length(xL);
shift = -min(xL)+log(n^(1/(n-1)));
Ay = [Ay;1];
Ax = [Ax;-ones(1,n)];
b = [b;log(exp(-shift)) + n*(shift)];