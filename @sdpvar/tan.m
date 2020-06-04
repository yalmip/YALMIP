function varargout = tan (varargin)
%TAN (overloaded)

switch class(varargin{1})

    case 'sdpvar'
        varargout{1} = InstantiateElementWise(mfilename,varargin{:});

    case 'char'

        operator = CreateBasicOperator('callback');        
        operator.bounds = @bounds;
        operator.derivative = @(x)(sec(x).^2);

        varargout{1} = [];
        varargout{2} = operator;
        varargout{3} = varargin{3};

    otherwise
        error(['SDPVAR/' upper(mfilename) ' called with weird argument']);
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