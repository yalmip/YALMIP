function varargout = SQRT(varargin)
%SQRT (overloaded)
%
% t = sqrt(x)
%
% The variable t can only be used in concavity preserving
% operations such as t>1, max t etc.
%
% When SQRT is used in a problem, the domain constraint
% (x>=0) is automatically added to the problem.
%
% See also SDPVAR, SDPVAR/GEOMEAN

switch class(varargin{1})
    
case 'double' % What is the numerical value of this argument (needed for displays etc)        
    X = sqrt(varargin{1});

case 'sdpvar' % Overloaded operator for SDPVAR objects. Pass on args and save them.
    X = varargin{1};
    [n,m] = size(X);
    if is(varargin{1},'real') & (n*m==1)
        varargout{1} = yalmip('addextendedvariable',mfilename,varargin{:});    
    else
        error('SQRT can only be applied to real scalars');
    end
    
case 'char' % YALMIP send 'model' when it wants the epigraph or hypograph
    if isequal(varargin{1},'graph')
        t = varargin{2}; % Second arg is the extended operator variable
        X = varargin{3}; % Third arg and above are the args user used when defining t.           
        varargout{1} = (cone([(X-1)/2;t],(X+1)/2));
        varargout{2} = struct('convexity','concave','monotonicity','increasing','definiteness','positive');
        varargout{3} = X;
    else           
    end    
otherwise
end
