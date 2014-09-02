function varargout = erfcx(varargin)
%ERFCX (overloaded)

switch class(varargin{1})

    case 'double'
        error('Overloaded SDPVAR/ERFCX CALLED WITH DOUBLE. Report error')

    case 'sdpvar'
        varargout{1} = InstantiateElementWise(mfilename,varargin{:});

    case 'char'

        operator = struct('convexity','none','monotonicity','decreasing','definiteness','positive','model','callback');
        operator.bounds = @bounds;

        varargout{1} = [];
        varargout{2} = operator;
        varargout{3} = varargin{3};

    otherwise
        error('SDPVAR/ERFCX called with CHAR argument?');
end

function [L,U] = bounds(xL,xU)
L = erfcx(xU);
U = erfcx(xL);