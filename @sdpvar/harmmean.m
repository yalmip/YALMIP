function varargout = harmmean(varargin)
%GEOMEAN (overloaded)
%
% t = HARMMEAN(X)
% For real vector X, returns length(X)/sum(X.^-1))

% See also SDPVAR, GEOMEAN

switch class(varargin{1})

    case 'sdpvar' % Overloaded operator for SDPVAR objects. Pass on args and save them.
        
        if nargin == 2 && isequal(varargin{2},2)
            varargout{1} = harmmean(varargin{1}')';
            return
        end
        
        X = varargin{1};
        [n,m] = size(X);
        if is(varargin{1},'hermitian') | min(n,m)==1
            varargout{1} = yalmip('define',mfilename,varargin{:});
        else
            % Create one variable for each column
            y = [];
            for i = 1:m
                index = (1+n*(i-1)):i*n;
                x = extsubsref(X,index);
                y = [y yalmip('define',mfilename,x)];
            end
            varargout{1} = y;
        end

    case 'char' % YALMIP send 'model' when it wants the epigraph or hypograph
        if isequal(varargin{1},'graph')
            t = varargin{2}; % Second arg is the extended operator variable
            X = varargin{3}; % Third arg and above are the args user used when defining t.

            n = length(X);
            z = sdpvar(n,1);          
            
            F = [X >= 0, sum(z) <= n*t];
            F = [F, cone([(X+z)';(X-z)';2*repmat(t,1,n)])];
            
            varargout{1} = F;
            varargout{2} = struct('convexity','concave','monotonicity','increasing','definiteness','positive');          
            varargout{3} = X;
        else

            varargout{1} = [];
            varargout{2} = [];
            varargout{3} = [];

        end
    otherwise
        error(['SDPVAR/' upper(mfilename) ' called with weird argument']);
end
