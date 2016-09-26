function varargout = interp1_internal(varargin)
%INTERP1_INTERNAL (overloaded)

switch class(varargin{1})

    case 'double'          
        varargout{1} = interp1(varargin{2},varargin{3},varargin{1},varargin{4});
            
    case 'char'

        operator = struct('convexity','none','monotonicity','none','definiteness','none','model','callback');
        operator.convexhull = [];
        operator.bounds = @bounds;

        varargout{1} = [min(varargin{4}) <= varargin{3} <= max(varargin{4})];
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
    L = min(yz);
    U = max(yz);    
end