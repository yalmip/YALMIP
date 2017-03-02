function varargout = interp1_internal(varargin)
%INTERP1_INTERNAL (overloaded)

switch class(varargin{1})

    case 'double'        
        varargout{1} = interp1(varargin{2},varargin{3},varargin{1},varargin{4});
      
    case 'char'
        
        t = varargin{2};
        x = varargin{3};
        xi = varargin{4};
        yi = varargin{5};
        method = varargin{6};
        
        if strcmpi(varargin{end},'milp')                        
            lambda = sdpvar(length(xi),1)            
            Model = [sos(lambda), x == lambda'*xi(:), t == lambda'*yi(:),lambda>=0, sum(lambda)==1];
            operator = struct('convexity','none','monotonicity','none','definiteness','none','model','integer');
        else
            Model = [min(xi) <= varargin{3} <= max(xi)];
            operator = struct('convexity','none','monotonicity','none','definiteness','none','model','callback');          
            operator.bounds = @bounds;
        end            
        varargout{1} = Model;
        varargout{2} = operator;
        varargout{3} = varargin{3};
        
    otherwise
        error('SDPVAR/INTERP1_INTERNAL called with strange argument');
end

function [L,U] = bounds(xL,xU,varargin)

xv = varargin{1};
yv = varargin{2};
if xL <= xv(1)
    index1 = 1;
else
    index1 = find(xL > xv);index1 = max(1,max(index1));
end
if xU >= xv(end)
    index2 = length(xv);
else
    index2 = find(xv < xU);index2 = min(length(xv),max(index2)+1);
end

if isequal(varargin{3},'linear')
    Z = yv(index1:index2);
    L = min(yv(index1:index2));
    U = max(yv(index1:index2));
else
    N = ceil((xv(index2)-xv(index1))/(mean(diff(xv))/100));
    z = linspace(xv(index1),xv(index2),N);
    yz = interp1(xv,yv,z,varargin{3});   
    % To account for finite grid, we add a precision dependent margin
    dz = (z(2)-z(1));
    L = min(yz)-dz;
    U = max(yz)+dz;    
end