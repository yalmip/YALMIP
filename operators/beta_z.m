function varargout=beta_z(varargin)

switch class(varargin{1})

    case 'double'
        varargout{1} = beta(varargin{1},varargin{2});
        
    case 'char'
        
        w = varargin{4};
        operator = CreateBasicOperator('convex','decreasing','positive','callback');
        operator.domain = [0 inf];        
        varargout{1} = [];
        varargout{2} = operator;
        varargout{3} = varargin{3};

    otherwise
        error(['SDPVAR/' upper(mfilename) ' called with weird argument']);
end