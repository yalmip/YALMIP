function varargout = sinh(varargin)

switch class(varargin{1})

    case 'sdpvar'
        varargout{1} = InstantiateElementWise(mfilename,varargin{:});
 
    case 'char'

        operator = CreateBasicOperator('increasing','callback');        
        operator.convexity = @convexity;
        operator.bounds = @bounds;
        operator.derivative = @(x)(cosh(x));

        varargout{1} = [];
        varargout{2} = operator;
        varargout{3} = varargin{3};

    otherwise
        error(['SDPVAR/' upper(mfilename) ' called with weird argument']);
end

function vexity = convexity(xL,xU)
if xL >= 0  
    vexity = 'convex';
elseif xU <= 0
    vexity = 'concave';
else
    vexity = 'none';
end