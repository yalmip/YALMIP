function varargout = pnorm(varargin)
%PNORM P-Norm of SDPVAR variable with convexity knowledge
%
% PNORM is used if your goal is to obtain a convex conic model model, since
% the function PNORM is implemented using conic epigrapraphs.
% (For p/q ==1,2,inf you should use the overloaded norm)
%
% t = pnorm(x,p/q), p/q >= 1
%
% By default pnorm is implemented using SOCP-represented cpower, which adds
% a large number of variables and constraints, in particular if p and q are
% large.
%
% An alternative is available if you have a power cone solver
%
% t = pnorm(x,p,'power')
%
% Note that this model is not based on a rational representation of the power

switch class(varargin{1})
    
    case {'double', 'gem', 'sgem'}
        varargout{1} = sum(abs(varargin{1}).^varargin{2}).^(1/varargin{2});
        
    case 'sdpvar' % Overloaded operator for SDPVAR objects. Pass on args and save them.
        X = varargin{1};
        [n,m] = size(X);
        if isreal(X) & min(n,m)==1
            if varargin{2}>=1
                varargout{1} = yalmip('define',mfilename,varargin{:});
            else
                error('PNORM only applicable for p>=1');
            end
        else
            error('PNORM can only be applied to real vectors.');
        end
        
    case 'char' % YALMIP send 'model' when it wants the epigraph or hypograph
        if isequal(varargin{1},'graph')
            t = varargin{2}; % Second arg is the extended operator variable
            X = varargin{3}; % Third arg and above are the args user used when defining t.
            p = varargin{4};
            if nargin > 4
                method = varargin{5};
            else
                method = 'socp';
            end
            
            switch method
                case 'power'
                    
                    n = length(X);
                    r = sdpvar(1,n);
                    F = [sum(r) == t, 
                         pcone([r;
                                repmat(t,1,n);
                                X(:)';
                                repmat(1/p,1,n)])];
                    
                case 'socp'
                    % sum(|xi|^(l/m))^(l/m) < t
                    % si>0
                    % -t^({l-m}/m)si^(m/l) < -xi
                    % -t^({l-m}/m)si^(m/l) <  xi
                    %sum(si) < t
                    % t>0                                                            
                    [l,m] = rat(p);
                    % l=4
                    % m = 1
                    % x<(t^(l-m)*s^m)^1/l
                    % -x<(t^(l-m)*s^m)^1/l
                    % x < t t t s
                    if 2^fix(log2(l))==l &  m == 1
                        s=sdpvar(length(X),1);
                        absX = sdpvar(length(X),1);
                        F = [-absX <= X <= absX];
                        F = [F,sum(s)<= t,s>=0];
                        for i = 1:length(X)
                            F = [F,detset(absX(i),[repmat(t,1,l-m) s(i)])];
                            F = [F,detset(-absX(i),[repmat(t,1,l-m) s(i)])];
                        end
                    else
                        % l = 7
                        % m = 2
                        % x < (t t t t t s s)^(1/7)
                        % x^7 < (t t t t t s s)
                        % x^8 < (t t t t t  s s x)
                        % x <  (t t t t t  s s x)^(1/8
                        s=sdpvar(length(X),1);
                        w = 2^(ceil(log2(l)));
                        absX = sdpvar(length(X),1);
                        F = [-absX <= X <= absX];
                        F = [F,sum(s)<= t,s>=0];
                        
                        for i = 1:length(X)
                            F = [F,detset(absX(i),[repmat(t,1,l-m) repmat(s(i),1,m) repmat(absX(i),1,w-l)])];
                            F = [F,detset(-absX(i),[repmat(t,1,l-m) repmat(s(i),1,m) repmat(absX(i),1,w-l)])];
                        end
                    end
                otherwise
                    error('Unknown method selected in pnorm (should be power or socp)');
            end
            
            varargout{1} = F;
            varargout{2} = struct('convexity','convex','monotonicity','none','definiteness','positive','model','graph');
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
