function varargout = atan2_internal(varargin)

switch class(varargin{1})

    case 'double' 
        
        varargout{1} = atan2(varargin{1}(1),varargin{1}(2));

    case 'char'
        
        operator = CreateBasicOperator('callback');
        
        varargout{1} = [];
        varargout{2} = operator;
        varargout{3} = varargin{3};
        
    otherwise
        error('SDPVAR/ATAN2 called with CHAR argument?');
end
