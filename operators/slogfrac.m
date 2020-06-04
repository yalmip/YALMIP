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
        if all(x==0)
            varargout{1} = log(2);% ?definition...
        elseif x(1)==0
            varargout{1} = log(1);
        else
            aux = 1 + (x(1)./(x(2)+sqrt(eps)));
            if aux <= 0
                varargout{1} = -15;
            else
                varargout{1} = log(1 + (x(1)./(x(2)+sqrt(eps))));
            end
            
        end
     
    case 'sdpvar'
        if min(size(varargin{1}))>1
            error('SLOGFRAC only defined for vector arguments');
        else
            varargout{1} = yalmip('define',mfilename,varargin{1});
        end

    case 'char'            
                                   
        operator = CreateBasicOperator('callback');
        operator.range = [-inf inf];
        operator.domain = [-inf inf];
        operator.bounds = @bounds;        
        operator.derivative = @(x) ([1./(x(1)+x(2)+sqrt(eps));1./(x(1)+x(2)+sqrt(eps))-(x(2)+sqrt(eps)).^-1]);
        
        varargout{1} = [];
        varargout{2} = operator;
        varargout{3} = varargin{3};

    otherwise
        error(['SDPVAR/' upper(mfilename) ' called with weird argument']);
end

function [L,U] = bounds(xL,xU)

% Derive bounds on x1/x2
z = [xL(1)./xL(2) xU(1)./xL(2) xL(1)./xU(2) xU(1)./xU(2)];
fL = min(z);
fU = max(z);

if fL <= -1
    L = -inf;
else
    L = log(1 + fL);
end
if fU <= -1
    U = -inf;
else
    U = log(1 + fU);
end


