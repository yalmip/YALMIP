function varargout = power_internal1(varargin)
%power_internal1
% Used for cases such as 2^x, and is treated as evaluation-based operators

switch class(varargin{1})

    case 'double'
        varargout{1} = varargin{2}.^varargin{1};

    case 'sdpvar'
        if isa(varargin{2},'sdpvar')
            x = varargin{2};
            y = varargin{1};
            varargout{1} = exp(y*log(x)); %x^y = exp(log(x^y))          
        else
            if length(varargin{1}) > 1 || size(varargin{2},1) ~= size(varargin{2},2)
                error('Inputs must be a scalar and a square matrix. To compute elementwise POWER, use POWER (.^) instead.');
            end
            x = varargin{2};
            y = varargin{1};
            if isa(x,'double') && x==1 && length(y)==1
                varargout{1} = 1;
            else                
                varargout{1} = InstantiateElementWise(mfilename,varargin{:});
            end
        end

    case 'char'
        
        X = varargin{3};
        Y = varargin{4};
        F=[];
        if Y>=1
            operator = CreateBasicOperator('increasing','positive','callback');
        elseif Y>=0
            operator = CreateBasicOperator('decreasing','positive','callback');
        else
            % Base is negative, so the power has to be an integer
            F = (integer(X));
            operator = CreateBasicOperator('decreasing','callback');
        end

        operator.bounds = @bounds_power;
        operator.convexhull = @convexhull_power;
        operator.derivative = @(x)derivative(x,Y);
        if Y >= 0
            operator.inverse = @(x,Y)inverse(x,Y);
        end

        varargout{1} = F;
        varargout{2} = operator;
        varargout{3} = [X(:);Y(:)];
    otherwise
        error([upper(mfilename) ' called with weird argument']);
end

% This should not be hidden here....
function [L,U] = bounds_power(xL,xU,base)
if base >= 1
    L = base^xL;
    U = base^xU;
elseif base>= 0
    L = base^xU;
    U = base^xL;
else
    disp('Not implemented yet. Report bug if you need this')
    error
end

function x = inverse(y,base)
if y <=0
    x = -inf;
else
    x = log(y)/log(base);
end

function df = derivative(x,base)
if length(base)~=length(x)
    base = base*ones(size(x));
end
f = base.^x;
df = log(base)*f;

function [Ax, Ay, b, K] = convexhull_power(xL,xU,base)
fL = base^xL;
fU = base^xU;
dfL = log(base)*fL;
dfU = log(base)*fU;
[Ax,Ay,b,K] = convexhullConvex(xL,xU,fL,fU,dfL,dfU);