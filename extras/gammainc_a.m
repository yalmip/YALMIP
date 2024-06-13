function varargout = gammainc_a(varargin)
%GAMMAINC_A

switch class(varargin{1})

    case 'double'       
        varargout{1} = gammainc(varargin{2},varargin{1});
  
    case 'char'

        operator = CreateBasicOperator('decreasing','positive','callback');            
        operator.range = [0 1];
        operator.domain = [0 inf];
        
        varargout{1} = [];
        varargout{2} = operator;
        varargout{3} = varargin{3};

    otherwise
        error(['SDPVAR/' upper(mfilename) ' called with weird argument']);
end
