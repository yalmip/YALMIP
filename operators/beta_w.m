function varargout=beta_w(varargin)

switch class(varargin{1})

    case 'double'
        varargout{1} = beta(varargin{2},varargin{1});
        
    case 'char'        
        operator = CreateBasicOperator('convex','decreasing','positive','callback');
        operator.domain = [0 inf];
        varargout{1} = [];
        varargout{2} = operator;
        varargout{3} = varargin{3};

    otherwise
        error(['SDPVAR/' upper(mfilename) ' called with weird argument']);
end