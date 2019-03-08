function varargout = atan2_internal(varargin)

switch class(varargin{1})

    case 'double' 
        varargout{1} = atan2(varargin{1}(1),varargin{1}(2));

    case 'char'
        X = varargin{3};
        varargout{1} = [];
        varargout{2} = struct('convexity','none','monotonicity','none','definiteness','none','model','callback');
        varargout{3} = [X(:)];
    otherwise
        error('SDPVAR/ATAN2 called with CHAR argument?');
end
