function varargout = gammaincinv_x(varargin)

switch class(varargin{1})

    case 'double'       
        varargout{1} = gammaincinv(varargin{1},varargin{2});
  
    case 'char'

        operator = CreateBasicOperator('increasing','positive','callback');        
        operator.range = [0 inf];
        operator.domain = [0 1];
        
        varargout{1} = [];
        varargout{2} = operator;
        varargout{3} = varargin{3};

    otherwise
        error(['SDPVAR/' upper(mfilename) ' called with weird argument']);
end
