function varargout = mvncdf(varargin)

switch class(varargin{1})

    case {'sdpvar','ndsdpvar'}
 
        if nargin > 1
            error('sdpvar/mvncdf currently only supports 1 argument, i.e. assumed zero mean and unit variance');
        end
        varargin{1} = reshape(varargin{1},[],1);
        varargout{1} = yalmip('define',mfilename,varargin{1});        

    case 'char'

        X = varargin{3};
        F = [];

        operator = CreateBasicOperator('positive','increasing','callback');  
        if length(X) == 1
            operator.convexity = @convexity; 
            operator.derivative = @(x)(1/2)*(1/sqrt(2))*exp(-(x/sqrt(2)).^2)*2/sqrt(pi);
            operator.inverse = @(x)(sqrt(2)*erfinv(2*x-1));
        else
            operator.bounds = @bounds;
        end
        
        operator.range = [0 1];        
        varargout{1} = F;
        varargout{2} = operator;
        varargout{3} = X;

    otherwise
        error(['SDPVAR/' upper(mfilename) ' called with weird argument']);
end

function [L,U] = bounds(xL,xU)
L = mvncdf(xL);
U = mvncdf(xU);

function vexity = convexity(xL,xU)
if xL >= 0  
    vexity = 'concave';
elseif xU <= 0
    vexity = 'convex';
else
    vexity = 'none';
end