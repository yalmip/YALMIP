function varargout = erfcinv(varargin)

switch class(varargin{1})

    case 'sdpvar'
        varargout{1} = InstantiateElementWise(mfilename,varargin{:});

    case 'char'
            
        operator = CreateBasicOperator('decreasing','callback');
        operator.convexity = @convexity;        
        operator.inverse = @(x)(erfc(x));
        operator.derivative = @(x)-1./(exp(-((erfcinv(x))).^2)*2/sqrt(pi));
        operator.domain = [-1 1];
        
        varargout{1} = [];
        varargout{2} = operator;
        varargout{3} = varargin{3};

    otherwise
        error(['SDPVAR/' upper(mfilename) ' called with weird argument']);
end

function vexity = convexity(xL,xU)
if xL >= 1 
    vexity = 'concave';
elseif xU <= 1
    vexity = 'convex';
else
    vexity = 'none';
end