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

switch class(varargin{1})

    case 'double'

        if nargin < 2
            M = repmat(1,size(varargin{1}));
        else
            M = varargin{2};
        end

        x = varargin{1};
        if isscalar(M) && ~isscalar(x)
            M = repmat(M,size(x));
        end
       
        y = x.^2;
        r = find(abs(x) > M);
        y(r) = M(r).*(2*abs(x(r)) - M(r));
        varargout{1} = sum(y);
        
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
        
        if ~isequal(size(M),size(X))
            M = repmat(M,size(X));
            if ~isequal(size(M),size(X))
                error('M must be scalar or same size as x')
            end
            varargin{2} = M;
        end
        
        y = [];
        for i = 1:size(X,2)
            varargin{1} = X(:,i);
            varargin{2} = M(:,i);
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

            % From Mangasarian & Musicant
            E = [-v <= X-u <= v, (u'*u + 2*sum(M.*v) <= t):'placeholder'];
            replacer.this = t;
            replacer.with = u'*u + 2*sum(M.*v);
            varargout{1} = E;
            varargout{2} = struct('convexity','convex','monotonicity','none','definiteness','positive','model','graph','replace',replacer);
            varargout{3} = X;

        else
            varargout{1} = [];
            varargout{2} = [];
            varargout{3} = [];
        end
    otherwise
end