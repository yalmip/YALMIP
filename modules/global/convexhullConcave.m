function [Ax,Ay,b,K] = convexhullConcave(varargin)

%function [Ax,Ay,b, K] = convexhullConcave(xL,xU,fL,fU,dfL,dfU)
% Two upper bounds from tangents
% y < f(xL) + (x-xL)*df(xL)
% y < f(xU) + (x-xL)*df(xU)
% lower bound from connecting extreme points
% y > f(xU)(x-xL)/(xU-xL) +  f(xL)(xU-x)/(xU-xL)
% can be written as
% Ax*x + Ay*y <= b

if rem(nargin,3)
    error('The convex hull generator assumes n triplets (x,f,df).')
end
m = nargin/3;
x = [varargin{(1:m)}]';
f = [varargin{(1:m)+m}]';
df = [varargin{(1:m)+2*m}]';

if df(1) == df(end)
    % This is a linear function
    % return y = f(1) + df*(x-x(1))
    Ax = df(1);
    Ay = -1;
    b = df(1)*x(1)-f(1);
    K.f = 1;
    K.l = 0;
    return
end

if all(diff(x)>=0) | all(diff(x(~isinf(x))))>=0 % Support [0 inf inf]
    Ay = -[-ones(m,1);1];
    b  = [f - x.*df; f(end)*x(1)/(x(end)-x(1)) -  f(1)*x(end)/(x(end)-x(1))];
    Ax  = [-df;f(end)/(x(end)-x(1)) - f(1)/(x(end)-x(1))];
    
    % Don't use ill-conditioned cuts
    if df(1)>1000
        Ax(1)=[];
        Ay(1) = [];
        b(1) = [];
    end
    if df(end)<-1000
        Ax(end-1)=[];
        Ay(end-1) = [];
        b(end-1) = [];
    end  
else
    Ax = [];
    Ay = [];
    b = [];    
end
j = find(any(isnan([Ax Ay b]),2));
if ~isempty(j)
    b(j)=[];
    Ax(j,:)=[];
    Ay(j,:)=[];
end
K.f = 0;
K.l = length(b);
if 0
  sdpvar x y;plot(Ax*x+Ay*y <= b);
end