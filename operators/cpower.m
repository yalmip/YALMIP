function varargout = cpower(varargin)
%CPOWER Power of SDPVAR variable with convexity knowledge
%
% CPOWER is recommended if your goal is to obtain
% a convex model, since the function CPOWER is implemented
% as a so called nonlinear operator. (For p/q ==2 you can
% however just as well use the overloaded power)
%
% t = cpower(x,p/q)
%
% For negative p/q, the operator is convex.
% For positive p/q with p>q, the operator is convex.
% For positive p/q with p<q, the operator is concave.
%
% A domain constraint x>0 is automatically added if
% p/q not is an even integer.
%
% Note, the complexity of generating the conic representation
% of these variables are O(2^L) where L typically is the
% smallest integer such that 2^L >= min(p,q)

switch class(varargin{1})

    case 'double'
        varargout{1} = power(varargin{1},varargin{2});
    case 'sdpvar' % Overloaded operator for SDPVAR objects. Pass on args and save them.
        X = varargin{1};
        if isreal(X)
            dim = size(X);
            X = reshape(X,prod(dim),1);
            y = [];
            for i = 1:prod(dim)
                y = [y;yalmip('define',mfilename,extsubsref(X,i),varargin{2})];
            end
            y = reshape(y,dim);
            varargout{1} = y;
        else
            error('CPOWER can only be applied to real vectors.');
        end

    case 'char' % YALMIP send 'model' when it wants the epigraph or hypograph
        if isequal(varargin{1},'graph')
            t = varargin{2}; % Second arg is the extended operator variable
            X = varargin{3}; % Third arg and above are the args user used when defining t.
            p = varargin{4};

            if p>0
                [p,q] = rat(abs(p));
                F = pospower(X,t,p,q);
                if p>q
                    convexity = 'convex';
                    monotonicity = 'increasing';
                else
                    convexity = 'concave';
                    monotonicity = 'decreasing';
                end
            else
                [p,q] = rat(abs(p));
                F = negpower(X,t,p,q);
                convexity = 'convex';
                monotonicity = 'decreasing';
            end

            varargout{1} = F;
            varargout{2} = struct('convexity',convexity,'monotonicity',monotonicity,'definiteness','positive','model','graph');
            varargout{3} = X;
        end
    otherwise
end

function F = pospower(x,t,p,q)
if p>q
    l = ceil(log2(abs(p)));
    r = 2^l-p;
    y = [ones(r,1)*x;ones(q,1)*t;ones(2^l-r-q,1)];
    F = detset(x,y);
else
    l = ceil(log2(abs(q)));
    y = [ones(p,1)*x;ones(2^l-q,1)*t;ones(q-p,1)];
    F = detset(t,y);
end

function F = negpower(x,t,p,q)
l = ceil(log2(abs(p+q)));
p = abs(p);
q = abs(q);
y = [ones(2^l-p-q,1);ones(p,1)*x;ones(q,1)*t];
F = detset(1,y);

