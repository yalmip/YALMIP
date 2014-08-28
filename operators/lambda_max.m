function varargout = lambda_max(varargin)
% LAMBDA_MAX Returns largest eigenvalue of Hermitian matrix.
%
% r = LAMBDA_MAX(X)
%
% See also LAMBDA_MIN

switch class(varargin{1})

    case 'double' % What is the numerical value of this argument (needed for displays etc)
        error(nargchk(1,1,nargin));
        X = varargin{1};
        [n,m] = size(X);
        if ~ishermitian(X)
            error('LAMBDA_MAX can only be applied on Hermitian matrices');
        else
            varargout{1} = max(real(eig(X)));
        end

    case 'sdpvar' % Overloaded operator for SDPVAR objects. Pass on args and save them.
        error(nargchk(1,1,nargin));
        X = varargin{1};
        [n,m] = size(X);
        if ~ishermitian(X)
            error('LAMBDA_MAX can only be Hermitian matrices');
        else
            varargout{1} = yalmip('define',mfilename,varargin{:});
        end

    case 'char' % YALMIP sends 'model' when it wants the epigraph or hypograph
        if isequal(varargin{1},'graph')
            t = varargin{2}; % Second arg is the extended operator variable
            X = varargin{3}; % Third arg and above are the args user used when defining t.
            varargout{1} = (t*eye(size(X,1)) >= X);
            varargout{2} = struct('convexity','convex','monotonicity','none','definiteness','none','model','graph');
            varargout{3} = X;
        else
            varargout{1} = [];
            varargout{2} = [];
            varargout{3} = [];
        end
    otherwise
end
