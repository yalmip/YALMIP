function varargout = cot(varargin)
%COT (overloaded)

switch class(varargin{1})

    case 'double'
        error('Overloaded SDPVAR/COT CALLED WITH DOUBLE. Report error')

    case 'sdpvar'
        varargout{1} = InstantiateElementWise(mfilename,varargin{:});

    case 'char'

        operator = struct('convexity','none','monotonicity','none','definiteness','none','model','callback');
        operator.derivative = @(x)(1./sin(x));

        varargout{1} = [];
        varargout{2} = operator;
        varargout{3} = varargin{3};

    otherwise
        error('SDPVAR/COS called with CHAR argument?');
end