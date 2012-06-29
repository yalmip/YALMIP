function varargout = invsathub (varargin)
%SIN (overloaded)

% Author Johan Löfberg
% $Id: invsathub.m,v 1.4 2008-02-18 15:42:45 joloef Exp $
switch class(varargin{1})

    case 'double'
        error('Overloaded SDPVAR/SIN CALLED WITH DOUBLE. Report error')

    case 'sdpvar'
        varargout{1} = InstantiateElementWise(mfilename,varargin{:});

    case 'char'

        operator = struct('convexity','none','monotonicity','none','definiteness','positive','model','callback');
        operator.bounds     = @bounds;
        operator.convexhull = @convexhull;
        %operator.derivative = @derivative;
        %operator.range = [-1 1];

        varargout{1} = [];
        varargout{2} = operator;
        varargout{3} = varargin{3};

    otherwise
        error('SDPVAR/SIN called with CHAR argument?');
end

function [L,U] = bounds(xL,xU,lambda)

fL = invsathub(xL,lambda);
fU = invsathub(xU,lambda);
U = max(fL,fU);
L = min(fL,fU);
if xL<0 & xU>0
    L = 0;
end

function [Ax, Ay, b, K] = convexhull(xL,xU,lambda)

K.l = 0;
K.f = 0;
fL = invsathub(xL,lambda);
fU = invsathub(xU,lambda);
if xU < -3*lambda % 1
    Ax = 0;Ay = 1;b = 3*lambda;K.f = 1;
elseif xL>=-lambda & xU <= 0 % 2
    Ax = 1;Ay = 1; b = 0;K.f = 1;
elseif xL>=0 & xU <= lambda % 3
    Ax = 1;Ay = -1; b = 0;K.f = 1;
elseif xU<=0 | xL>=0 %4
    dfL = derivative(xL,lambda);
    dfU = derivative(xU,lambda);
    [Ax,Ay,b,K] = convexhullConcave(xL,xU,fL,fU,dfL,dfU);
elseif xL<0 & xL>=-lambda & xU>0 & xU<= lambda % 5
    [Ax,Ay,b,K] = convexhullConvex(xL,xU,fL,fU,-lambda,lambda);
elseif 0%xL<=-lambda & xL>=-3*lambda
    z = [xL 0 xU];
    fz = [fL 0 fU];
    k1 = max((fz(2:end)-fz(1))./(z(2:end)-xL))+1e-12;
    k2 = min((fz(2:end)-fz(1))./(z(2:end)-xL))-1e-12;
    k3 = min((fz(1:end-1)-fz(end))./(z(1:end-1)-xU))+1e-12;
    k4 = max((fz(1:end-1)-fz(end))./(z(1:end-1)-xU))-1e-12;
    Ax = [-k1;k2;-k3;k4];
    Ay = [1;-1;1;-1];
    b =  [k1*(-z(1)) + fz(1);-(k2*(-z(1)) + fz(1));k3*(-z(end)) + fz(end);-(k4*(-z(end)) + fz(end))];
    K.l = length(b);
elseif xL<-3*lambda & xU>3*lambda
     z = [xL 0 xU];
    fz = [fL 0 fU];
    k1 = max((fz(2:end)-fz(1))./(z(2:end)-xL))+1e-12;
    k2 = min((fz(2:end)-fz(1))./(z(2:end)-xL))-1e-12;
    k3 = min((fz(1:end-1)-fz(end))./(z(1:end-1)-xU))+1e-12;
    k4 = max((fz(1:end-1)-fz(end))./(z(1:end-1)-xU))-1e-12;
    Ax = [-k1;k2;-k3;k4];
    Ay = [1;-1;1;-1];
    b =  [k1*(-z(1)) + fz(1);-(k2*(-z(1)) + fz(1));k3*(-z(end)) + fz(end);-(k4*(-z(end)) + fz(end))];
    K.l = length(b);
  
%       clf
%          x = linspace(xL,xU,1000);
%          plot(polytope([Ax Ay],b));   hold on
%          plot(x,invsathub(x,lambda))
%          1

else
    
    z = [linspace(xL,xU,100)];
    fz = [invsathub(z,lambda)];
    if xU>0 & xL<0
        z = [0 z];
        fz = [0 fz];
    end

    [minval,minpos] = min(fz);
    [maxval,maxpos] = max(fz);
    xtestmin = linspace(z(max([1 minpos-5])),z(min([100 minpos+5])),100);
    xtestmax = linspace(z(max([1 maxpos-5])),z(min([100 maxpos+5])),100);
    fz1 = invsathub(xtestmin,lambda);
    fz2 = invsathub(xtestmax,lambda);
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
    K.l = length(b);

end

% 
%   clf
%          x = linspace(xL,xU,1000);
%          plot(polytope([Ax Ay],b));   hold on
%          plot(x,invsathub(x,lambda))
%          1
         
function df=derivative(x,lambda)
if nargin==1
    lambda=0.5;
end
df = 0;
if (-3*lambda < x) & (x < -lambda)
    df = -0.25*(2*x+6*lambda);
elseif (-lambda < x) & (x < 0)
    df = -lambda;
elseif (x>0) & (x<lambda)
    df = lambda;
elseif (x>lambda) & (x<3*lambda)
    df = -0.25*(2*x-6*lambda);
end