function varargout = asec(varargin)
%ACOT (overloaded)

% Author Johan Löfberg
% $Id: asec.m,v 1.4 2007-08-02 18:16:26 joloef Exp $
switch class(varargin{1})

    case 'double'
        error('Overloaded SDPVAR/ACOT CALLED WITH DOUBLE. Report error')

    case 'sdpvar'
        varargout{1} = InstantiateElementWise(mfilename,varargin{:});

    case 'char'

        operator = struct('convexity','none','monotonicity','none','definiteness','none','model','callback');
        operator.convexhull = [];
        operator.bounds = @bounds;
        operator.derivative = @(x)(1./(x.*(x.^2-1).^0.5));

        varargout{1} = [];
        varargout{2} = operator;
        varargout{3} = varargin{3};

    otherwise
        error('SDPVAR/ACOT called with CHAR argument?');
end

function [L,U] = bounds(xL,xU)
if xL<=-1 & xU >=1
    L = -inf;
    U = inf;
else
    L = asec(xU);
    U = asec(xL);
end