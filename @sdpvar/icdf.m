function varargout = icdf(varargin)
%ICDF (overloaded)

switch class(varargin{2})    
    case 'sdpvar'
        temp = {varargin{2},varargin{1},varargin{3:end}};
        varargout{1} = InstantiateElementWise('icdf_internal',temp{:});        
    otherwise        
end