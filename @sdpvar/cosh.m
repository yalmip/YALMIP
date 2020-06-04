function varargout = cosh(varargin)
%COSH (overloaded)

switch class(varargin{1})

    case 'sdpvar'
        varargout{1} = InstantiateElementWise(mfilename,varargin{:});

    case 'char'

        operator = CreateBasicOperator('convex','positive','callback');        
        operator.bounds = @bounds;
        operator.derivative = @(x)(sinh(x));

        varargout{1} = [];
        varargout{2} = operator;
        varargout{3} = varargin{3};

    otherwise
        error(['SDPVAR/' upper(mfilename) ' called with weird argument']);
end

function [L,U] = bounds(xL,xU)
if xL<0 & xU>0
    L = 0;
    U = max([cosh(xL) cosh(xU)]);
elseif xL<0
    L = cosh(xU);
    U = cosh(xL);
else
    L = cosh(xL);
    U = cosh(xU);
end
