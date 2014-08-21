function varargout = gammaln(varargin)
%GAMMALN (overloaded)

switch class(varargin{1})

    case 'sdpvar'
        varargout{1} = InstantiateElementWise(mfilename,varargin{:});

    case 'char'

        operator = struct('convexity','convex','monotonicity','none','definiteness','none','model','callback');
       
        operator.domain = [0 inf];      
        operator.derivative =@(x)psi(0,x);                
        
        varargout{1} = [];
        varargout{2} = operator;
        varargout{3} = varargin{3};

    otherwise
        error('SDPVAR/GAMMALN called with CHAR argument?');
end