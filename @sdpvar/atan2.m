function varargout = atan2(varargin)

switch class(varargin{1})

    case 'sdpvar' 
        varargout{1} = yalmip('define','atan2_internal',[varargin{:}]);        

    case 'char'

        X = varargin{3};
        
        operator = CreateBasicOperator('callback');
        operator.range = [-pi pi];
        
        varargout{1} = [];
        varargout{2} = operator;
        varargout{3} = [X(:)];
    otherwise
        error('SDPVAR/atan2 called with CHAR argument?');
end
