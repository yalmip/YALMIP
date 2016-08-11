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

        X = varargin{3};
       
        operator = struct('convexity','convex','monotonicity','none','definiteness','none','model','callback');        
        operator.bounds = @bounds;
        operator.derivative = @(x) exp(x)./(sum(exp(x)));
        operator.convexhull = @convexhull;
                
        varargout{1} = [];
        varargout{2} = operator;
        varargout{3} = X;

    otherwise
        error('SDPVAR/LOG called with CHAR argument?');
end


function [L, U] = bounds(xL,xU)

L = log(sum(exp(xL)));
U = log(sum(exp(xU)));

function [Ax, Ay, b] = convexhull(xL,xU)
Ax = [];
Ay = [];
b = [];