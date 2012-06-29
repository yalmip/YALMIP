function varargout = acosh(varargin)
%ACOSH (overloaded)

% Author Johan Löfberg
% $Id: acosh.m,v 1.4 2007-08-02 18:16:26 joloef Exp $
switch class(varargin{1})

    case 'double'
        error('Overloaded SDPVAR/ASIN CALLED WITH DOUBLE. Report error')

    case 'sdpvar'
        varargout{1} = InstantiateElementWise(mfilename,varargin{:});

    case 'char'

        operator = struct('convexity','none','monotonicity','none','definiteness','positive','model','callback');
        operator.convexhull = [];
        operator.bounds = @bounds;

        varargout{1} = [];
        varargout{2} = operator;
        varargout{3} = varargin{3};

    otherwise
        error('SDPVAR/ACOSH called with CHAR argument?');
end

function [L,U] = bounds(xL,xU)
if xU<=-1
    L = asinh(xU);
    U = asinh(xL);
elseif xL>=1
    L = asinh(xL);
    U = asinh(xU);
else
    L = -inf;
    U = inf;
end

