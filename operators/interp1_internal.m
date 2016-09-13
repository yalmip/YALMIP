function varargout = interp1_internal(varargin)
%INTERP1_INTERNAL (overloaded)

switch class(varargin{1})

    case 'double'   
        varargout{1} = interp1(varargin{2},varargin{3},varargin{1});
            
    case 'char'

        operator = struct('convexity','none','monotonicity','none','definiteness','none','model','callback');
        operator.convexhull = [];
        operator.bounds = @bounds;

        varargout{1} = [];
        varargout{2} = operator;
        varargout{3} = varargin{3};
         
    otherwise
        error('SDPVAR/INTERP1_INTERNAL called with strange argument');
end

function [L,U] = bounds(xL,xU,varargin)

xv = varargin{1};
yv = varargin{2};
[xL xU]
if xL <= xv(1)
    index1 = 1;
else
    index1 = find(L > xv);index1 = max(1,min(index1)-1);
end

if xU >= xv(end)
    index2 = length(xv);
else
index2 = find(xv < xU);index2 = min(length(xv),max(index1)+1);
end
L = min(yv(index1:index2));
U = max(yv(index1:index2));