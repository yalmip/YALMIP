function varargout = gammainc_x(varargin)
%GAMMAINC_X

switch class(varargin{1})

    case 'double'       
        varargout{1} = gammainc(varargin{1},varargin{2});
  
    case 'char'

        operator = CreateBasicOperator('increasing','positive','callback');        
        operator.range = [0 1];
        operator.domain = [0 inf];
        
        varargout{1} = [];
        varargout{2} = operator;
        varargout{3} = varargin{3};

    otherwise
        error('SDPVAR/GAMMAINC_X called with strange argument?');
end
