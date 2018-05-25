function varargout = interp1_nonlinear(varargin)
%INTERP1_NONLINEAR (overloaded)

switch class(varargin{1})
    
    case 'double'
    % YALMIP shuffles around argument and places y-values (to be optimized)
    % in first argument  
    
    if isempty(varargin{2})
        xifi = varargin{1};
        xi = xifi(1:length(xifi)/2);
        fi = xifi(length(xifi)/2+1:end);
        varargout{1} = interp1(xi,fi,varargin{3},varargin{4});         
    else
        varargout{1} = interp1(varargin{1},varargin{2},varargin{3},varargin{4});
    end
        
    case 'char' 
        if isempty(varargin{4})
            xifi = varargin{3};
            xi = xifi(1:length(xifi)/2);
            F = diff(xi) >= 0;
        else
            F = [];
        end
        operator = struct('convexity','none','monotonicity','none','definiteness','none','model','callback');                
        varargout{1} = F;
        varargout{2} = operator;
        varargout{3} = varargin{3};
        
    otherwise
        error('SDPVAR/INTERP1_NONLINEAR called with strange argument');
end