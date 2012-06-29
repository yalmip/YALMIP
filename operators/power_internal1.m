function varargout = power_internal(varargin)
%power_internal1
% Used for cases such as 2^x, and is treated as evaluation-based operators

% Author Johan Löfberg
% $Id: power_internal1.m,v 1.7 2007-08-07 11:16:18 joloef Exp $
switch class(varargin{1})

    case 'double'
        varargout{1} = varargin{2}.^varargin{1};

    case 'sdpvar'
        if isa(varargin{2},'sdpvar')
            error('x^y currently not supported for SDPVAR x and SDPVAR y')
        else
            varargout{1} = InstantiateElementWise(mfilename,varargin{:});
        end

    case 'char'
        
        X = varargin{3};
        Y = varargin{4};
        F=[];
        if Y>=1
            operator = struct('convexity','none','monotonicity','increasing','definiteness','positve','model','callback');
        elseif Y>=0
            operator = struct('convexity','none','monotonicity','decreasing','definiteness','positive','model','callback');
        else
            % Base is negative, so the power has to be an integer
            F = set(integer(X));
            operator = struct('convexity','none','monotonicity','decreasing','definiteness','none','model','callback');
        end

        operator.bounds = @bounds_power;
        operator.convexhull = @convexhull_power;

        varargout{1} = F;
        varargout{2} = operator;
        varargout{3} = [X(:);Y(:)];
    otherwise
        error('SDPVAR/power_internal1 called with CHAR argument?');
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

function [Ax, Ay, b] = convexhull_power(xL,xU,base)
fL = base^xL;
fU = base^xU;
dfL = log(base)*fL;
dfU = log(base)*fU;
[Ax,Ay,b] = convexhullConvex(xL,xU,fL,fU,dfL,dfU);