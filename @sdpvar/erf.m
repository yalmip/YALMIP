function varargout = erf(varargin)
%ERF (overloaded)

switch class(varargin{1})

    case 'sdpvar'
        varargout{1} = InstantiateElementWise(mfilename,varargin{:});

    case 'char'

        operator = CreateBasicOperator('increasing','callback');
        operator.convexity = @convexity;        
        operator.range = [-1 1];
        operator.derivative =@(x)exp(-x.^2)*2/sqrt(pi);
        operator.inverse = @(x)(erfinv(x));
              
        varargout{1} = [];
        varargout{2} = operator;
        varargout{3} = varargin{3};

    otherwise
        error(['SDPVAR/' upper(mfilename) ' called with weird argument']);
end

function vexity = convexity(xL,xU)
if xL >= 0  
    vexity = 'concave';
elseif xU <= 0
    vexity = 'convex';
else
    vexity = 'none';
end

