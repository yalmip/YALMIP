function varargout = sumk(varargin)
% SUMK  Returns sum of k largest (eigen-)values.
%
% s = SUMK(X,k)
%
% For a vector X, SUMK returns the sum of the k largest elements.
%
% For a symmetric matrix X, SUMK returns the sum of the k largest eigen-values.
%
% See also SUMABSK

% ***************************************************
% This file defines a nonlinear operator for YALMIP
%
% It can take three different inputs
% For double inputs, it returns standard double values
% For sdpvar inputs, it genreates a an internal variable
% When first input is 'model' it generates the epigraph
% and a structure descring properties of the operator
% ***************************************************
switch class(varargin{1})

    case 'double' % What is the numerical value of this argument (needed for displays etc)
        if nargin == 1
            error('sumk needs two arguments');
        else
            X = varargin{1};
            [n,m] = size(X);
            if (min(n,m) > 1 & ~issymmetric(X))
                error('sumk can only be applied on vectors and symmetric matrices');
            else
                k = min(length(X),varargin{2});
                if min(n,m)==1
                    sorted = sort(X);
                else
                    sorted = sort(eig(X));
                end
                varargout{1} = sum(sorted(max(1,end-k+1):end));
            end
        end

    case 'sdpvar' % Overloaded operator for SDPVAR objects. Pass on args and save them.
        X = varargin{1};
        [n,m] = size(X);
        if (min(n,m) > 1 & ~issymmetric(X))
            error('sumk can only be applied on vectors and symmetric matrices');
        else
            if nargin<2
                error('sumk needs two arguments');
            else
                varargout{1} = yalmip('define',mfilename,varargin{:});
            end
        end

    case 'char' % YALMIP send 'model' when it wants the epigraph or hypograph
        if isequal(varargin{1},'graph')
            t = varargin{2}; % Second arg is the extended operator variable
            X = varargin{3}; % Third arg and above are the args user used when defining t.
            k = min(varargin{4},length(X));
            [n,m] = size(X);
            Z = sdpvar(n,m);
            s = sdpvar(1,1);
            if min(n,m)==1
                varargout{1} = (t-k*s-sum(Z) >= 0) + (Z >= 0) + (Z-X+s >= 0);
                varargout{2} = struct('convexity','convex','monotonicity','increasing','definiteness','none','model','graph');
            else
                varargout{1} = (t-k*s-trace(Z) >= 0) + (Z >= 0) + (Z-X+s*eye(n) >= 0);
                varargout{2} = struct('convexity','convex','monotonicity','none','definiteness','none','model','graph');
            end
            varargout{3} = X;
        else
        end
    otherwise
end