function varargout = sec(varargin)
%SEC (overloaded)

switch class(varargin{1})

    case 'double'
        error('Overloaded SDPVAR/SEC CALLED WITH DOUBLE. Report error')

    case 'sdpvar'
        varargout{1} = InstantiateElementWise(mfilename,varargin{:});

    case 'char'

        operator = struct('convexity','none','monotonicity','none','definiteness','none','model','callback');
        operator.derivative = @(x)(tan(x).*sec(x));

        varargout{1} = [];
        varargout{2} = operator;
        varargout{3} = varargin{3};

    otherwise
        error('SDPVAR/SEC called with CHAR argument?');
end