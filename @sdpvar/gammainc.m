function varargout = gammainc(varargin)

if nargin ~= 2
    error('Not enough input arguments in gammainc.');
end

if isa(varargin{1},'sdpvar') && isa(varargin{2},'sdpvar')
    error('Only one argument in gammainc(X,A) can be an SDPVAR');
end

switch class(varargin{1})
    case 'double'
        
        if varargin{1} < 0
            error('X must be real and non-negative in gammainc(X,A)')
        end
        varargout{1} = InstantiateElementWise('gammainc_a',varargin{2:-1:1});
        
    case 'sdpvar'
        
        if varargin{2}<0
            error('A must be real and non-negative gammainc(X,A)');
        end
        varargout{1} = InstantiateElementWise('gammainc_x',varargin{:});
        
    otherwise
        error(['SDPVAR/' upper(mfilename) ' called with weird argument']);
end