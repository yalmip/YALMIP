function varargout = sec(varargin)

switch class(varargin{1})

    case 'sdpvar'
        varargout{1} = InstantiateElementWise(mfilename,varargin{:});

    case 'char'

        operator = CreateBasicOperator('callback');
        operator.derivative = @(x)(tan(x).*sec(x));
        operator.convexity = @convexity;        
        operator.inversebounds = @inversebounds;

        varargout{1} = [];
        varargout{2} = operator;
        varargout{3} = varargin{3};

    otherwise
        error(['SDPVAR/' upper(mfilename) ' called with weird argument']);
end

function vexity = convexity(xL,xU)
if xL >= -pi/2 && xU <= pi/2
    vexity = 'convex';
else
    vexity = 'none';
end

function [xLi,xUi] = inversebounds(fL,fU,xL,xU)

if xL <= 0 && xU >= 0 && xL >=-pi/2 && xU <= pi/2
    % Convex region
    a = asec(fU);
    xUi = a;
    xLi = -a;
elseif xL>=0 && xU <= pi/2
    % Increasing region
    xLi = asec(fL);
    xUi = asec(fU);
elseif xL >=-pi/2 && xU <= 0
    % Decreasing region
    xLi = -asec(fU);
    xUi = -asec(fL);
else
    xLi = xL;
    xUi = xU;
end