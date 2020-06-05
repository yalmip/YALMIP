function varargout=beta(varargin)

switch class(varargin{1})

    case 'sdpvar'
        
        if ~isa(varargin{2},'double')
            error('W is not allowed to be an SDPVAR object')
        end
        varargout{1} = InstantiateElementWise(mfilename,varargin{:});

    case 'char'
        
        operator = CreateBasicOperator('convex','decreasing','positive','callback');
        operator.domain = [0 inf];
        varargout{1} = [];
        varargout{2} = operator;
        varargout{3} = varargin{3};

    otherwise
        error(['SDPVAR/' upper(mfilename) ' called with weird argument']);
end