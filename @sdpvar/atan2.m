function varargout = atan2(varargin)

switch class(varargin{1})

    case 'sdpvar' 
        varargout{1} = yalmip('define','atan2_internal',[varargin{:}]);        
   
    otherwise
        error('SDPVAR/atan2 called with CHAR argument?');
end