function varargout = norm_nuclear(varargin)
% NORM_NUCLEAR Returns sum of singular values
%
% s = SUMK(X,k)
%
% For a vector X, NORM_NUCLEAR returns the sum of singular values
%
% For a matrix X, NORM_NUCLEAR returns the sum of absolute values.
%
% See also SUMABSK

switch class(varargin{1})
    
    case 'double' % What is the numerical value of this argument (needed for displays etc)
        varargout{1} = sum(svd(varargin{1}));
        
    case 'sdpvar'
        varargout{1} = yalmip('define',mfilename,varargin{:});
        
    case 'char' % YALMIP sends 'model' when it wants the epigraph or hypograph
        if isequal(varargin{1},'graph')
            t = varargin{2};
            X = varargin{3};
            [n,m] = size(X);
            if is(X,'real')
                U = sdpvar(m);
                V = sdpvar(n);
            else
                U = sdpvar(m,m,'hermitian','complex');
                V = sdpvar(n,n,'hermitian','complex');
            end
            F = [trace(U)+trace(V) <= 2*t, [U X';X V]>=0];
            varargout{1} = F;
            varargout{2} = struct('convexity','convex','monotonicity','none','definiteness','positive','model','graph');
            varargout{3} = X;
        else
            varargout{1} = [];
            varargout{2} = [];
            varargout{3} = [];
        end
    otherwise
end