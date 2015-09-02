function varargout = kullbackleibler(varargin)
% KULLBACKLEIBLER
%
% y = KULLBACKLEIBLER(x,y)
%
% Computes/declares the convex Kullback-Leibler divergence sum(x.*log(x./y))
% Alternatively -sum(x.*log(y/x)), i.e., negated perspectives of log(y)
%
% See also ENTROPY, CROSSENTROPY


switch class(varargin{1})

    case 'double'    
        if nargin == 1
            z = varargin{1};
            z = reshape(z,[],2);
            x = z(:,1);
            y = z(:,2);
        else
            x = varargin{1}(:);
            y = varargin{2}(:);
        end
        l = log(x./y);
        l(x<=0) = 0;
        l = real(l);
        varargout{1} = sum(x.*l);       

    case {'sdpvar','ndsdpvar'}

        varargin{1} = reshape(varargin{1},[],1);
        varargin{2} = reshape(varargin{2},[],1);    
        
        if length(varargin{1})~=length(varargin{2})
            if length(varargin{1})==1
                varargin{1} = repmat(varargin{1},length(varargin{2}),1);
            elseif  length(varargin{2})==1
                varargin{2} = repmat(varargin{2},length(varargin{1}),1);
            else
                error('Dimension mismatch in kullbackleibler')
            end
        end
        
        varargout{1} = yalmip('define',mfilename,[varargin{1};varargin{2}]);
        
    case 'char'

        X = varargin{3};
        F = [X >= 0];

        operator = struct('convexity','convex','monotonicity','none','definiteness','none','model','callback');
        operator.range = [-inf inf];
        operator.domain = [0 inf];       
        operator.derivative = @derivative;

        varargout{1} = F;
        varargout{2} = operator;
        varargout{3} = X;

    otherwise
        error('SDPVAR/KULLBACKLEIBLER called with CHAR argument?');
end

function df = derivative(x)
z = abs(reshape(x,[],2));
x = z(:,1);
y = z(:,2);
% Use KL = -Entropy + Cross Entropy
df = [1+log(x);0*y] + [-log(y);-x./y];


