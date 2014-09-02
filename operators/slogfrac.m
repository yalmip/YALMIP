function varargout = slogfrac(varargin)
%logfrac
%
% y = SLOGFRAC(x)
%
% Computes/declares log(1 + x(1)/x(2))

switch class(varargin{1})
    case 'double'
        x = varargin{1};        
        % Safe version with defined negative values (helps fmincon when
        % outside feasible region)  
        if all(x)==0
            varargout{1} = log(2);% ?definition...
        elseif x(1)==0
            varargout{1} = log(1);
        else
            varargout{1} = log(1 + (x(1)./(x(2)+sqrt(eps))));
        end
     
    case 'sdpvar'
        if min(size(varargin{1}))>1
            error('SLOGFRAC only defined for vector arguments');
        else
            varargout{1} = yalmip('define',mfilename,varargin{1});
        end

    case 'char'            
        X = varargin{3};
        Y = X(2);
        X = X(1);             
        F = [];%X > sqrt(eps), Y > sqrt(eps)];
        operator = struct('convexity','none','monotonicity','none','definiteness','none','model','callback');
        operator.range = [0 inf];
        operator.domain = [1e-6 inf];
        operator.bounds = @bounds;        
        operator.derivative = @(x) ([1./(x(1)+x(2)+sqrt(eps));1./(x(1)+x(2)+sqrt(eps))-(x(2)+sqrt(eps)).^-1]);
        varargout{1} = F;
        varargout{2} = operator;
        varargout{3} = varargin{3};

    otherwise
        error('SDPVAR/SLOGFRAC called with CHAR argument?');
end

