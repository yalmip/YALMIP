function varargout = log2(varargin)
%log2 (overloaded)

varargout{1} = log(varargin{1})/log(2);
return

% Author Johan Löfberg
% $Id: log2.m,v 1.7 2007-08-02 18:16:26 joloef Exp $
switch class(varargin{1})

    case 'double' 
        error('Overloaded SDPVAR/LOG2 CALLED WITH DOUBLE. Report error')

    case 'sdpvar'
        varargout{1} = check_for_special_cases(varargin{:});
        % Nope, then just define this logarithm
        if isempty(varargout{1})
            varargout{1} = InstantiateElementWise(mfilename,varargin{:});
        end
        
    case 'char'      

        X = varargin{3};
        F = [X >= 1e-8];

        operator = struct('convexity','concave','monotonicity','increasing','definiteness','none','model','callback');
        operator.convexhull = @convexhull;
        operator.bounds = @bounds;
        operator.derivative = @derivative;

        varargout{1} = F;
        varargout{2} = operator;
        varargout{3} = X;

    otherwise
        error('SDPVAR/LOG2 called with CHAR argument?');
end

function df = derivative(x)
df = (1./x)/log(2);

function [L,U] = bounds(xL,xU)
if xL < 0
    % The variable is not bounded enough yet
    L = -inf;
elseif xL==0
    L = -inf;
else
    L = log2(xL);
end
if xU < 0
    % This is an infeasible problem
    L = inf;
    U = -inf;
else
    U = log2(xU);
end

function [Ax, Ay, b] = convexhull(xL,xU)
fL = log2(xL);
fU = log2(xU);
dfL = (1/(xL))/log(2);
dfU = (1/(xU))/log(2);
[Ax,Ay,b] = convexhullConcave(xL,xU,fL,fU,dfL,dfU);

function f = check_for_special_cases(x)
f = [];
% Check for log2(1+x)
base = getbase(x);
if all(base(:,1)==1)
    f = slog2(x-1);
    return;
end