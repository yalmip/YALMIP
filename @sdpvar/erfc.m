function varargout = erfc(varargin)
%ERFC (overloaded)

% Author Johan Löfberg
% $Id: erfc.m,v 1.6 2008-02-28 21:54:49 joloef Exp $
switch class(varargin{1})

    case 'double'
        error('Overloaded SDPVAR/ERFC CALLED WITH DOUBLE. Report error')

    case 'sdpvar'
        varargout{1} = InstantiateElementWise(mfilename,varargin{:});

    case 'char'

        operator = struct('convexity','none','monotonicity','decreasing','definiteness','positive','model','callback');
        operator.bounds = @bounds;
        operator.range = [-1 1];
        operator.derivative =@(x)-exp(-x.^2)*2/sqrt(pi);
        varargout{1} = [];
        varargout{2} = operator;
        varargout{3} = varargin{3};

    otherwise
        error('SDPVAR/ERF called with CHAR argument?');
end

function [L,U] = bounds(xL,xU)
L = erfc(xL);
U = erfc(xU);