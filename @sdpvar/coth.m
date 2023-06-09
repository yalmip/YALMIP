function varargout = coth(varargin)

switch class(varargin{1})

    case 'sdpvar'
        varargout{1} = InstantiateElementWise(mfilename,varargin{:});

    case 'char'

        operator = CreateBasicOperator('callback');
        operator.monotonicity = @monotonicity;                                      
        operator.derivative = @(x)(-(1-coth(x)).^2); 
        operator.singularity = [0 -inf inf];
        operator.inflection = [-inf -1 0 1];
               
        varargout{1} = [];
        varargout{2} = operator;
        varargout{3} = varargin{3};

    otherwise
        error(['SDPVAR/' upper(mfilename) ' called with weird argument']);
end

function mono = monotonicity(xL,xU)
if xL >= 0 || xU <= 0
    mono = 'decreasing';
else
    mono = 'none';
end