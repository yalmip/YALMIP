function varargout = erfcinv(varargin)
%ERFCINV (overloaded)

switch class(varargin{1})

    case 'double'
        error('Overloaded SDPVAR/ERFCINV CALLED WITH DOUBLE. Report error')

    case 'sdpvar'
        varargout{1} = InstantiateElementWise(mfilename,varargin{:});

    case 'char'

        operator = struct('convexity','none','monotonicity','decreasing','definiteness','none','model','callback');
        operator.bounds = @bounds;

        varargout{1} = [];
        varargout{2} = operator;
        varargout{3} = varargin{3};

    otherwise
        error('SDPVAR/ERFCINV called with CHAR argument?');
end

function [L,U] = bounds(xL,xU)
L = erfcinv(xU);
U = erfcinv(xL);