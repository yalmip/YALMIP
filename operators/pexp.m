function varargout = pexp(varargin)
%PEXP
%
% y = PEXP(x)
%
% Computes perspective exp, x(1)*exp(x(2)/x(1)) on x>0
%
% Implemented as evalutation based nonlinear operator. Hence, the convexity
% of this function is exploited to perform convexity analysis and rigorous
% modelling.

switch class(varargin{1})
    
    case 'double'
        
        if ~isequal(prod(size(varargin{1})),2)
            error('PEXP only defined for 2x1 arguments');
        end
        x = varargin{1};
        
        varargout{1} = x(1)*exp(x(2)/x(1));


    case 'sdpvar'

        if ~isequal(prod(size(varargin{1})),2)
            error('PLOG only defined for 2x1 arguments');
        else
            varargout{1} = yalmip('define',mfilename,varargin{1});
        end

    case 'char'

        X = varargin{3};

        operator = struct('convexity','convex','monotonicity','none','definiteness','positive','model','callback');
        operator.range = [0 inf];
        operator.domain = [0 inf];    
        operator.derivative = @derivative;

        varargout{1} = [];
        varargout{2} = operator;
        varargout{3} = X;

    otherwise
        error('SDPVAR/PEXP called with CHAR argument?');
end

function dp = derivative(x)
z = x(2)/x(1);
dp = [exp(z)-z*exp(z);exp(z)];


