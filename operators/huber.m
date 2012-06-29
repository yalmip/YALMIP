function varargout = huber(varargin)
% HUBER  Returns the Hüber function (convex operator)
%
% p = HUBER(x,M)
%
% Returns sum_i h(x_i,M) where h(x) is the scalar Hüber function
%  h(x) = x^2         if |x|<=M
%         M(2*|x|-M)  otherwise
%
% The Hüber functions is convex and non-monotonic.

% Author Johan Löfberg
% $Id: huber.m,v 1.6 2009-05-29 08:05:12 joloef Exp $

switch class(varargin{1})

    case 'double'

        if nargin < 2
            M = 1;
        else
            M = varargin{2};
        end

        x = varargin{1};
        y = x.^2;
        r = find(abs(x) > M);
        y(r) = M.*(2*abs(x(r)) - M);
        varargout{1} = y;
        
    case 'sdpvar' % Pass on args and save them.

        if nargin < 2
            M = 1;
            varargin{end+1} = 1;
        else
            M = varargin{2};
        end

        X = varargin{1};
        [n,m] = size(X);
        if min(n,m) == 1
            X = X(:);
        end
        y = [];
        for i = size(X,2);
            varargin{1} = X(:,i);
            y = [y yalmip('define',mfilename,varargin{:})];
        end
        varargout{1} = y;

    case 'char' % YALMIP send 'model' when it wants the epigraph or hypograph
        if isequal(varargin{1},'graph')
            t = varargin{2};
            X = varargin{3};
            M = varargin{4};

            u = sdpvar(length(X),1);
            v = sdpvar(length(X),1);

            % From Boyd & Vandenberghe
%            E = [v>=0, 0 <= u <= M, -u-v <= X <= u + v, u'*u + 2*M*sum(v) <= t];
            % From Mangasarian & Musicant
            E = [-v <= X-u <= v, u'*u + 2*M*sum(v) <= t];
            
            varargout{1} = E;
            varargout{2} = struct('convexity','convex','monotonicity','none','definiteness','positive','model','graph');
            varargout{3} = X;

        else
            varargout{1} = [];
            varargout{2} = [];
            varargout{3} = [];
        end
    otherwise
end