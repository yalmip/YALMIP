function varargout = atan2_internal(varargin)

switch class(varargin{1})

    case 'double' 
        
        varargout{1} = atan2(varargin{1}(1),varargin{1}(2));

    case 'char'
        
        operator = CreateBasicOperator('callback');
        operator.derivative = @(x)derivative(x);
        operator.range = [-pi pi];
        
        varargout{1} = [];
        varargout{2} = operator;
        varargout{3} = varargin{3};
        
    otherwise
        error([upper(mfilename) ' called with weird argument']);
end

function df = derivative(x)

% z = x(1)/x(2), f = atan(z)
% f' = z'*1/(1+z^2)
if ~any(x)
    % Undefined really..
    df = [0 ; 0];
elseif x(2)==0
    df = [0 ; -1/x(1)];
else
    z = x(1)/x(2);
    dzdx1 = 1/x(2);
    dzdx2 = -x(1)/x(2)^2;
    df = [dzdx1 ; dzdx2]/(1+z^2);
end