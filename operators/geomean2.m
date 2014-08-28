function varargout = geomean2(varargin)
%GEOMEAN2 Nonlinear operator in YALMIP
%
% t = GEOMEAN2(X)
%
% For Hermitian matrix X, returns det(X)^(1/(2^ceil(log2(length(X)))))
%
% For real vector X, returns prod(X)^(1/(2^ceil(log2(length(X)))))
%
% This concave function is monotonically growing in det(P)
% for P>0, so  it can be used for maximizing det(P), 
% or to add lower bound constraints on the determinant.
%
% When GEOMEAN2 is used in a problem, the domain constraint
% X>=0 is automatically added to the problem.
%
% Note that the function is the geometric mean of
% the elements (or eigenvalues) if the dimension of 
% X is a power of 2, hence the name GEOMEAN2.
%
% See also SDPVAR, SDPVAR/GEOMEAN, SUMK, SUMABSK

switch class(varargin{1})
    
case 'double' % What is the numerical value of this argument (needed for displays etc)        
    X = varargin{1};
    [n,m] = size(X);
    if min(n,m)==1
        varargout{1} = prod(X)^(1/(2^ceil(log2(length(X)))));
    else
        if norm(X-X')<n*eps
            varargout{1} = (real(det(X)))^(1/(2^ceil(log2(length(X)))));
        else
            error('GEOMEAN2 can only be applied to real vectors and Hermitian matrices');
        end
    end
    
case 'sdpvar' % Overloaded operator for SDPVAR objects. Pass on args and save them.
    X = varargin{1};
    [n,m] = size(X);
    if is(varargin{1},'hermitian') | min(n,m)==1
        varargout{1} = yalmip('define',mfilename,varargin{:});    
    else
        error('GEOMEAN2 can only be applied to real vectors and Hermitian matrices');
    end
    
case 'char' % YALMIP send 'model' when it wants the epigraph or hypograph
    if isequal(varargin{1},'graph')
        t = varargin{2}; % Second arg is the extended operator variable
        X = varargin{3}; % Third arg and above are the args user used when defining t.           
        varargout{1} = detset(t,X);        
        if issymmetric(X)
            varargout{2} = struct('convexity','concave','monotonicity','none','definiteness','positive','model','graph');
        else
            varargout{2} = struct('convexity','concave','monotonicity','increasing','definiteness','positive','model','graph');
        end
        varargout{3} = X;
    else           
    end    
otherwise
end
