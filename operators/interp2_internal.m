function varargout = interp2_internal(varargin)
%INTERP1_INTERNAL (overloaded)

switch class(varargin{1})

    case 'double'          
        varargout{1} = interp2(varargin{2},varargin{3},varargin{4},varargin{1}(1),varargin{1}(2),varargin{5});
            
    case 'char'        
        L = [min(varargin{4}(:));min(varargin{5}(:))];
        U = [max(varargin{4}(:));max(varargin{5}(:))];
        operator = struct('convexity','none','monotonicity','none','definiteness','none','model','callback');      
        operator.bounds = @bounds;
        varargout{1} = [L <= varargin{3} <= U];
        varargout{2} = operator;
        varargout{3} = varargin{3};         
    otherwise
        error('SDPVAR/INTERP1_INTERNAL called with strange argument');
end

function [L,U] = bounds(xL,xU,varargin)

xv = varargin{1};
yv = varargin{2};
zv = varargin{3};
i = 1;
lowindex = zeros(2,1);
highindex = zeros(2,1);
if xL(1) <= xv(1)
    lowindex(i) = 1;
else
    temp = find(xL(1) > xv(1,:));lowindex(i) = max(1,max(temp));
end
if xU(1) >= xv(end)
    highindex(i) = size(xv,2);
else
    temp = find(xv(1,:) < xU(1));highindex(i) = min(size(xv,2),max(temp)+1);
end
i = 2;
if xL(2) <= yv(1)
    lowindex(i) = 1;
else
    temp = find(xL(2) > yv(:,1));lowindex(i) = max(1,max(temp));
end
if xU(2) >= yv(end)
    highindex(i) = size(yv,1);
else
    temp = find(yv(:,1) < xU(2));highindex(i) = min(size(yv,1),max(temp)+1);
end

if isequal(varargin{4},'linear')
    Z = zv(lowindex(1):highindex(1),lowindex(2):highindex(2));
    L = min(Z(:));
    U = max(Z(:));
else   
    dx = (xv(1,highindex(1))-xv(1,lowindex(1)))/(size(xv,2)*25);
    dy = (yv(highindex(2),1)-yv(lowindex(2),1))/(size(yv,1)*25);    
    [X,Y] = meshgrid(xv(1,lowindex(1)):dx:xv(1,highindex(1)),yv(lowindex(2),1):dy:yv(highindex(2),1));
    Z = interp2(xv,yv,zv,X,Y,varargin{4});
    % To account for finite grid, we add a precision dependent margin
    L = min(Z(:))-dx-dy;
    U = max(Z(:))+dx+dy;
end