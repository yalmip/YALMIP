function varargout = erf(varargin)
%ERF (overloaded)

% Author Johan Löfberg
% $Id: erf.m,v 1.12 2008-05-06 16:37:50 joloef Exp $
switch class(varargin{1})

    case 'double'
        error('Overloaded SDPVAR/ERF CALLED WITH DOUBLE. Report error')

    case 'sdpvar'
        varargout{1} = InstantiateElementWise(mfilename,varargin{:});

    case 'char'

        operator = struct('convexity','none','monotonicity','increasing','definiteness','none','model','callback');
        operator.bounds = @bounds;
        operator.range = [-1 1];
        operator.derivative =@(x)exp(-x.^2)*2/sqrt(pi);
        operator.convexhull = @convexhull;
        
        varargout{1} = [];
        varargout{2} = operator;
        varargout{3} = varargin{3};

    otherwise
        error('SDPVAR/ERF called with CHAR argument?');
end

function [L,U] = bounds(xL,xU)
L = erf(xL);
U = erf(xU);


function [Ax, Ay, b, K] = convexhull(xL,xU)
K = [];
if xU <= 0    
    xM = (xL+xU)/2;
    fL = erf(xL);
    fM = erf(xM);
    fU = erf(xU);
    dfL = exp(-xL.^2)*2/sqrt(pi);
    dfM = exp(-xM.^2)*2/sqrt(pi);
    dfU = exp(-xU.^2)*2/sqrt(pi);
    
    [Ax,Ay,b] = convexhullConvex(xL,xM,xU,fL,fM,fU,dfL,dfM,dfU);
       
elseif xL >= 0
    xM = (xL+xU)/2;
    fL = erf(xL);
    fM = erf(xM);
    fU = erf(xU);
    dfL = exp(-xL.^2)*2/sqrt(pi);
    dfM = exp(-xM.^2)*2/sqrt(pi);
    dfU = exp(-xU.^2)*2/sqrt(pi);    
    [Ax,Ay,b] = convexhullConcave(xL,xM,xU,fL,fM,fU,dfL,dfM,dfU);
else
    z = linspace(xL,xU,1000);
    fz = erf(z);
    [minval,minpos] = min(fz);
    [maxval,maxpos] = max(fz);
    xtestmin = linspace(z(max([1 minpos-5])),z(min([100 minpos+5])),100);
    xtestmax = linspace(z(max([1 maxpos-5])),z(min([100 maxpos+5])),100);

    fz1 = erf(xtestmin);
    fz2 = erf(xtestmax);
    z = [z(:);xtestmin(:);xtestmax(:)];
    fz = [fz(:);fz1(:);fz2(:)];
    [z,sorter] = sort(z);
    fz = fz(sorter);
    [z,ii,jj]=unique(z);
    fz = fz(ii);

    k1 = max((fz(2:end)-fz(1))./(z(2:end)-xL))+1e-12;
    k2 = min((fz(2:end)-fz(1))./(z(2:end)-xL))-1e-12;
    k3 = min((fz(1:end-1)-fz(end))./(z(1:end-1)-xU))+1e-12;
    k4 = max((fz(1:end-1)-fz(end))./(z(1:end-1)-xU))-1e-12;
    Ax = [-k1;k2;-k3;k4];
    Ay = [1;-1;1;-1];
    b =  [k1*(-z(1)) + fz(1);-(k2*(-z(1)) + fz(1));k3*(-z(end)) + fz(end);-(k4*(-z(end)) + fz(end))];
end
   



