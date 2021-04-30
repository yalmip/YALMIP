function varargout = sumabsk(varargin)
% SUMK  Returns sum of k largest absolute (eigen-)values.
%
% s = SUMABSK(X,k,w)
%
% For a vector X, SUMABSK returns the sum of the k largest absolute elements.
%
% For a symmetric matrix X, SUMABSK returns the sum of the k largest
% absolute eigen-values. 
%
% A third argument can be used to generalize to weighted sum.
% The weights have to be non-negative and non-increasing.
%
% See also SUMABSK

if nargin == 1
    error('sumk needs at least two arguments');
end

switch class(varargin{1})
    
    case 'double' % What is the numerical value of this argument (needed for displays etc)
        
        X = varargin{1};
        [n,m] = size(X);
        if (min(n,m) > 1 & ~issymmetric(X))
            error('sumk can only be applied on vectors and symmetric matrices');
        else
            k = min(length(X),varargin{2});
            w = ones(k,1);
            if nargin == 3
                w = varargin{3};
                if length(w)==1
                    w = ones(k,1)*w;               
                end
            end
            if min(n,m)==1
                sorted = sort(X,'descend');
            else
                sorted = sort(eig(X),'descend');
            end
            sorted = sorted(:);
            w = w(:);
            varargout{1} = sum(sorted(1:k).*w(1:k));
        end
        
    case 'sdpvar' % Overloaded operator for SDPVAR objects. Pass on args and save them.
        X = varargin{1};
        [n,m] = size(X);
        if nargin < 3
            w = 1;
        else
            w = varargin{3};
            if length(w)>1
                if any(diff(w)>0) || any(w<0)
                    error('The weights have to be non-negative non-increasing')
                end
            end
        end
        if (min(n,m) > 1 & ~issymmetric(X))
            error('sumk can only be applied on vectors and symmetric matrices');
        else
            varargout{1} = yalmip('define',mfilename,varargin{:});
        end
        
    case 'char' % YALMIP send 'model' when it wants the epigraph or hypograph
        if isequal(varargin{1},'graph')
            t = varargin{2}; % Second arg is the extended operator variable
            X = varargin{3}; % Third arg and above are the args user used when defining t.
            k = min(varargin{4},length(X));
            if nargin == 5
                w = varargin{5};
            else
                w = 1;
            end                
            [Model,Properties] = sumabsk_generator(X,k,t,w);
            varargout{1} = Model;
            varargout{2} = Properties;
            varargout{3} = X;
        else
        end
    otherwise
end
