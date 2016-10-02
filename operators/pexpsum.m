function varargout = pexpsum(varargin)
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
        x = varargin{1};
        x = reshape(x,2,[]);
        varargout{1} = sum(x(1,:).*exp(x(2,:)./x(1,:)));


    case 'sdpvar'

        if ~(isequal(prod(size(varargin{1})),2) || size(varargin{1},1)==2)
            error('PEXP only defined for 2xN arguments');
        else
            x = reshape(varargin{1},[],1);            
            varargout{1} = yalmip('define',mfilename,varargin{1});                         
        end

    case 'char'

        X = varargin{3};

        operator = struct('convexity','convex','monotonicity','none','definiteness','positive','model','callback');
        operator.range = [0 inf];
        operator.domain = [0 inf];    
        operator.derivative = @derivative;
        operator.convexhull = @convexhull;
        operator.bounds = @bounds;
        
        varargout{1} = [];
        varargout{2} = operator;
        varargout{3} = X;

    otherwise
        error('SDPVAR/PEXP called with CHAR argument?');
end

function dp = derivative(x)
x = reshape(x,2,[]);
z = x(2,:)./x(1,:);
z(x(1,:)==0)=0;
dp = [exp(z)-z.*exp(z);exp(z)];
dp = dp(:);
