function varargout = acos(varargin)
%ACOS (overloaded)

switch class(varargin{1})

    case 'sdpvar'
        % We will compute value in an internal file in order to extend
        % the function beyond -1 and 1 (since fmincon makes calls outside
        % bounds...)
        varargout{1} = InstantiateElementWise('acos_internal',varargin{:});

    otherwise
        error('SDPVAR/ACOS called with CHAR argument?');
end
