function varargout = sdpfun(varargin)
%SDPFUN Gateway to general (elementwise) functions on SDPVAR variables (overloaded)

% Author Johan Löfberg
% $Id: sdpfun.m,v 1.12 2007-08-02 19:17:36 joloef Exp $

switch class(varargin{1})

    case {'struct','double'} % What is the numerical value of this argument (needed for displays etc)
        if isequal(varargin{end}(1),'@')
            fun = eval(varargin{end});
            varargin{end} = fun;
        end
        varargout{1} = feval(varargin{end},varargin{1:end-1});

    case 'sdpvar' % Overloaded operator for SDPVAR objects. Pass on args and save them.
        % try to figure out size of expected output (many->1 or many->many
        if isequal(varargin{end}(1),'@')
            fun = eval(varargin{end});
            varargin{end} = fun;
        end
        varargin_copy = varargin;
        varargin_copy{1} = zeros(size(varargin_copy{1},1),size(varargin_copy{1},2));
        output = feval(varargin_copy{end},varargin_copy{1:end-1});
        if prod(size(output)) == 1
            % MANY -> 1
            varargout{1} = yalmip('define',mfilename,varargin{:});
        else
            % MANY -> MANY
            y = [];
            [n,m] = size(varargin{1});
            varargin{1} = varargin{1}(:);
            for i = 1:length(varargin{1})
                varargin_copy = varargin;
                varargin_copy{1} = varargin_copy{1}(i);
                y = [y;yalmip('define',mfilename,varargin_copy{:})];
            end
            varargout{1} = reshape(y,n,m);
        end

    case 'char' % YALMIP sends 'model' when it wants the epigraph or hypograph
      
        operator = struct('convexity','none','monotonicity','none','definiteness','none','model','callback');
        operator.bounds = @bounds;
        operator.convexhull = @convexhull;

        varargout{1} = [];
        varargout{2} = operator;
        varargout{3} = varargin{3};

    otherwise
        error('SDPVAR/SDPFUN called with CHAR argument?');
end

function [Ax,Ay,b] = convexhull(xL,xU,varargin)
if ~(isinf(xL) | isinf(xU))
    z = linspace(xL,xU,100);
    fz = feval(varargin{end},z,varargin{1:end-1});
    % create 4 bounding planes
    % f(z) < k1*(x-XL) + f(xL)
    % f(z) > k2*(x-XL) + f(xL)
    % f(z) < k3*(x-XU) + f(xU)
    % f(z) > k4*(x-XU) + f(xU)
    % f(z) < k5*(x-XM) + f(xM)
    % f(z) > k6*(x-XM) + f(xM)
 
    k1 = max((fz(2:end)-fz(1))./(z(2:end)-xL))+1e-12;
    k2 = min((fz(2:end)-fz(1))./(z(2:end)-xL))-1e-12;
    k3 = min((fz(1:end-1)-fz(end))./(z(1:end-1)-xU))+1e-12;
    k4 = max((fz(1:end-1)-fz(end))./(z(1:end-1)-xU))-1e-12;
    Ax = [-k1;k2;-k3;k4];
    Ay = [1;-1;1;-1];
    b =  [k1*(-z(1)) + fz(1);-(k2*(-z(1)) + fz(1));k3*(-z(end)) + fz(end);-(k4*(-z(end)) + fz(end))];
else
    Ax = [];
    Ay = [];
    b = [];
end

function [L,U] = bounds(xL,xU,varargin)

if ~isinf(xL) & ~isinf(xU)
    xtest = linspace(xL,xU,100);
    values = feval(varargin{end},xtest,varargin{1:end-1});
    [minval,minpos] = min(values);
    [maxval,maxpos] = max(values);
    % Concetrate search around the extreme regions
    xtestmin = linspace(xtest(max([1 minpos-5])),xtest(min([100 minpos+5])),100);
    xtestmax = linspace(xtest(max([1 maxpos-5])),xtest(min([100 maxpos+5])),100);
    values1 = feval(varargin{end},xtestmin,varargin{1:end-1});
    values2 = feval(varargin{end},xtestmax,varargin{1:end-1});
    L = min([values1 values2]);
    U = max([values1 values2]);
else
    L = -inf;
    U = inf;
end
