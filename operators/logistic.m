function varargout = logistic(varargin)
% LOGISTIC  Returns logistic function 1./(1+exp(-x))
%
% y = LOGISTIC(x)
%
% For a real vector x, LOGISTIC returns (1+exp(-x)).^-1

switch class(varargin{1})

    case 'double'
        x = varargin{1};
        varargout{1} = 1./(1+exp(-x));

    case 'sdpvar'
        varargout{1} = InstantiateElementWise(mfilename,varargin{:});

    case 'char'
        
        operator = CreateBasicOperator('increasing','positive','callback');        
        operator.derivative = @(x)logistic(x).*(1-logistic(x));
        operator.inverse    = @(x)(log(x)-log(1-x));
        operator.inflection = [0 -1];
        operator.range = [0 1];

        varargout{1} = [];
        varargout{2} = operator;
        varargout{3} = varargin{3};

    otherwise
        error([upper(mfilename) ' called with weird argument']);
end