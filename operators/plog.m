function varargout = plog(varargin)
%PLOG
%
% y = PLOG(x)
%
% Computes concave perspective log, x(1)*log(x(2)/x(1)) on x>0
%
% Implemented as evalutation based nonlinear operator. Hence, the concavity
% of this function is exploited to perform convexity analysis and rigorous
% modelling.

switch class(varargin{1})
    
    case 'double'
        
        if ~isequal(prod(size(varargin{1})),2)
            error('PLOG only defined for 2x1 arguments');
        end
        x = varargin{1};
        % Safe version with defined negative values (helps fmincon when
        % outside feasible region)

        if isequal(x(1),[0])
            varargout{1} = 0;
        else
            varargout{1} = x(1)*log(x(2)/x(1));
        end

    case 'sdpvar'

        if ~isequal(prod(size(varargin{1})),2)
            error('PLOG only defined for 2x1 arguments');
        else
            varargout{1} = yalmip('define',mfilename,varargin{1});
        end

    case 'char'

        X = varargin{3};
      
        operator = struct('convexity','concave','monotonicity','none','definiteness','none','model','callback');
        operator.range = [-inf inf];
        operator.domain = [0 inf];
      %  operator.bounds = @bounds;
      %  operator.convexhull = @convexhull;
        operator.derivative = @derivative;

        varargout{1} = [];
        varargout{2} = operator;
        varargout{3} = X;

    otherwise
        error('SDPVAR/PLOG called with CHAR argument?');
end

function dp = derivative(x)
z = x(2)/x(1);
dp = [log(z)-1;1./z];

function [L,U] = bounds(xL,xU)
xU(isinf(xU)) = 1e12;
x1 = xL(1)*log(xL(1)/xL(2));
x2 = xU(1)*log(xU(1)/xL(2));
x3 = xL(1)*log(xL(1)/xU(2));
x4 = xU(1)*log(xU(1)/xU(2));
U = max([x1 x2 x3 x4]);
if (exp(-1)*xU(2) > xL(1)) & (exp(-1)*xU(2) < xU(1))
    L = -exp(-1)*xU(2);
else
    L = min([x1 x2 x3 x4]);
end

function [Ax,Ay,b] = convexhull(L,U);%xL,xU)
% 
% Ax = [];
% Ay = [];
% b = [];
% return

% L = [1 2];
% U = [3 6];

p1 = [L(1);L(2)];
p2 = [U(1);L(2)];
p3 = [L(1);U(2)];
p4 = [U(1);U(2)];
pm = [(L(1) + U(1))/2;(L(2) + U(2))/2];

z1 = L(1)*log(L(1)/L(2));
z2 = U(1)*log(U(1)/L(2));
z3 = L(1)*log(L(1)/U(2));
z4 = U(1)*log(U(1)/U(2));
zm = pm(1)*log(pm(1)/pm(2));

g1 = [log(L(1)/L(2)) + 1;-L(1)/L(2)];
g2 = [log(U(1)/L(2)) + 1;-U(1)/L(2)];
g3 = [log(L(1)/U(2)) + 1;-L(1)/U(2)];
g4 = [log(U(1)/U(2)) + 1;-U(1)/U(2)];
gm = [log(pm(1)/pm(2)) + 1;-pm(1)/pm(2)];
%sdpvar x y z
%C = [max([z1 z2 z3 z4]) > z > z1 + g1'*([x;y] - p1),z > z2 + g2'*([x;y] - p2), z > z3 + g3'*([x;y] - p3),z > z4 + g4'*([x;y] - p4),[x y]>U, L < [x y]]
Ax = [0 0;g1';g2';g3';g4';gm'];%;eye(2);-eye(2)];
Ay = [1;-1;-1;-1;-1;-1];%;0;0;0;0];
b = [max([z1 z2 z3 z4]);g1'*p1-z1;g2'*p2-z2;g3'*p3-z3;g4'*p4-z4;gm'*pm-zm];%;U(:);-L(:)];

xL = L(1);
yL = L(2);
xU = U(1);
yU = U(2);

p1 = [xL;yL;z1];
p2 = [xU;yL;z2];
p3 = [xL;yU;z3];
p4 = [xU;yU;z4];
U = max([z1 z2 z3 z4]);
H1 = null([p2-p1 p3-p1]')';
H2 = null([p3-p2 p4-p2]')';
A = [H1;H2];
b = [b;H1*p1;H2*p2];
Ax = [Ax;A(:,1:2)];
Ay = [Ay;A(:,3)];


