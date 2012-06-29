function varargout = absexp(varargin)
%LOG (overloaded)

% Author Johan Löfberg
% $Id: absexp.m,v 1.1 2007-08-01 15:03:13 joloef Exp $

switch class(varargin{1})

    case 'double'
        varargout{1} = abs(exp(varargin{1}) - 1);

    case 'sdpvar'
        varargout{1} = InstantiateBuiltInScalar(mfilename,varargin{:});

    case 'char'
        t = varargin{2};
        X = varargin{3};

        F = SetupEvaluationVariable(varargin{:});

        operator = struct('convexity','none','monotonicity','none','definiteness','none','model','callback');
   %     operator.convexhull = @convexhull;
        operator.bounds = @bounds;

        varargout{1} = F;
        varargout{2} = operator;
        varargout{3} = X;

    otherwise
        error('SDPVAR/LOG called with CHAR argument?');
end

function [L,U] = bounds(xL,xU)
if xL <= 0 & xU>=0
    % The variable is not bounded enough yet
    L = 0;
    U = max([abs(exp(xL)-1) abs(exp(xU)-1)]);
else
    U = max([abs(exp(xL)-1) abs(exp(xU)-1)]);
    L = min([abs(exp(xL)-1) abs(exp(xU)-1)]);
end

function [Ax, Ay, b] = convexhull(xL,xU)
Ax = [];
Ay = [];
b = [];
