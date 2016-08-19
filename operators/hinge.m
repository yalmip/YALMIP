function varargout = hinge(varargin)
%HINGE Models convex operator max(0,x^p) for integer p>=1
%
% t = hinge(x,p)


switch class(varargin{1})

    case 'double'
        if nargin == 1
            varargout{1} = max(0,varargin{1});
        else
            varargout{1} = max(0,varargin{1}.^varargin{2});
        end
        
    case 'sdpvar' % Overloaded operator for SDPVAR objects. Pass on args and save them.
        X = varargin{1};
        if nargin == 1
            p = 2;
        else
            p = varargin{2};
        end
        if p <= 1
            error('HINGE max(0,x^p) must have p>=1')
        elseif ~(p==ceil(p))
            error('HINGE max(0,x^p) must have integer p')
        elseif p == 1
            varargout{1} = max(0,X);
        else
            varargout{1} = yalmip('define',mfilename,X,p);
        end
        
    case 'char' % YALMIP send 'model' when it wants the epigraph or hypograph
        if isequal(varargin{1},'graph')
            t = varargin{2}; % Second arg is the extended operator variable
            X = varargin{3}; % Third arg and above are the args user used when defining t.
            p = varargin{4};

            convexity = 'convex';
            if ~even(p)                    
                    monotonicity = 'increasing';
            else                    
                    monotonicity = 'none';
            end
            
            e = sdpvar(1);
            F = [e >= X, pospower(e,t,p)];                
            
            varargout{1} = F;
            varargout{2} = struct('convexity',convexity,'monotonicity',monotonicity,'definiteness','positive','model','graph');
            varargout{3} = X;
        end
    otherwise
end

function F = pospower(x,t,p)
q = 1;
l = ceil(log2(abs(p)));
r = 2^l-p;
y = [ones(r,1)*x;ones(q,1)*t;ones(2^l-r-q,1)];
F = detset(x,y);

