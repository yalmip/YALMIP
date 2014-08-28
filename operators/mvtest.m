function varargout=mvtest(varargin)
%EIG (overloaded)

switch class(varargin{1})

    case 'double'
        varargout{1} = sum(varargin{1}.^2,2);

    case 'sdpvar' % Overloaded operator for SDPVAR objects. Pass on args and save them.

        x = varargin{1};
        y = yalmip('definemulti',mfilename,x,[size(x,1) 1]);
        varargout{1} = y;

    case 'char'

        X = varargin{3};

        operator = struct('convexity','convex','monotonicity','none','definiteness','none','model','callback');

        varargout{1} = [];
        varargout{2} = operator;
        varargout{3} = X;

    otherwise
        error('Strange type on first argument in SDPVAR/SORT');
end
