function varargout = tan (varargin)
%TAN (overloaded)

% Author Johan Löfberg
% $Id: tan.m,v 1.10 2007-08-02 18:16:26 joloef Exp $

switch class(varargin{1})

    case 'double'
        error('Overloaded SDPVAR/TAN CALLED WITH DOUBLE. Report error')

    case 'sdpvar'
        varargout{1} = InstantiateElementWise(mfilename,varargin{:});

    case 'char'

        operator = struct('convexity','none','monotonicity','none','definiteness','none','model','callback');
        operator.convexhull = [];
        operator.bounds = @bounds;
        operator.derivative = @(x)(-1./cos(x));

        varargout{1} = [];
        varargout{2} = operator;
        varargout{3} = varargin{3};

    otherwise
        error('SDPVAR/TAN called with CHAR argument?');
end


function [L,U] = bounds(xL,xU)
n1 = fix((xL+pi/2)/(pi));
n2 = fix((xU+pi/2)/(pi));
if n1==n2
    L = tan(xL);
    U = tan(xU);
else
    L = -inf;
    U = inf;
end